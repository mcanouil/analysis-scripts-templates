#' qc_sample_sheet_meqtl
#' @import data.table
qc_sample_sheet_meqtl <- function(phenotype, methylation, exclusion, relatedness, ethnicity) {
  methy_sample_sheet <- data.table::fread(file = methylation)

  relatedness <- data.table::fread(file = relatedness)

  related_pairs <- relatedness[
    i = sub("_.*", "", IID1) != sub("_.*", "", IID2)
  ]

  best_related_samples <- sapply(
    X = unique(lapply(
      X = related_pairs[["IID1"]],
      FUN = function(iid) {
        related_pairs[IID1 %in% iid | IID2 %in% iid, sort(unique(c(IID1, IID2, iid)))]
      }
    )),
    FUN = function(x) {
      data.table::melt.data.table(
        data = related_pairs[IID1 %in% x | IID2 %in% x, .SD, .SDcols = paste0(rep(c("IID", "F_MISS"), each = 2), 1:2)],
        measure.vars = patterns("^F_MISS", "^IID"),
        value.name = c("F_MISS", "IID")
      )[
        which.min(F_MISS),
        unique(IID)
      ]
    }
  )

  bad_duplicated_samples <- relatedness[
    i = sub("_.*", "", IID1) == sub("_.*", "", IID2)
  ][
    j = data.table::melt.data.table(.SD, measure.vars = paste0("F_MISS", 1:2)),
    .SDcols = c(paste0("IID", 1:2), paste0("F_MISS", 1:2))
  ][
    j = list(
      IID = sub("_.*", "", IID1),
      best_iid = data.table::fifelse(
        test = variable[which.min(value)] == "F_MISS1",
        yes = IID1,
        no = IID2
      )
    ),
    by = c(paste0("IID", 1:2))
  ]

  exclusion <- data.table::fread(file = exclusion)[
    i = Sex_Discrepancy == 1 | Sample_Call_Rate == 1 | Sex_Missing == 1 | Heterozygosity_Check == 1,
    j = Status := "Exclude"
  ][
    i = IID %in%
      bad_duplicated_samples[j = unlist(.SD, use.names = FALSE), .SDcols = patterns("IID[0-9]+")] &
        Status == "Check",
    j = Status := NA_character_
  ]

  dt <- phenotype[j = `:=`(
    "Sample_Name" = IID,
    "vcf_id" = IID
  )][# Add related variable
    j = is_related := vcf_id %in% relatedness[
      i = sub("_.*", "", IID1) != sub("_.*", "", IID2),
      j = unique(unlist(.SD)),
      .SDcols = paste0("IID", 1:2)
    ]
  ][
    j = is_best_related := vcf_id %in% best_related_samples
  ][# keep best duplicates from genotyping
    i = bad_duplicated_samples,
    j = vcf_id := best_iid,
    on = "IID"
  ]

  dt_snp <- merge(
    x = merge(x = dt, y = exclusion[j = list(vcf_id = IID, Status)], by = "vcf_id", all.x = TRUE),
    y = data.table::fread(file = ethnicity)[j = -c("cohort")], # Add genetics PCs
    by.x = "vcf_id",
    by.y = "iid",
    all.x = TRUE
  )[
    i = (is_related),
    j = Status := "Related"
  ][
    i = (is_related) & (is_best_related),
    j = Status := "Related (Best)"
  ][
    i = is.na(PC01),
    j = `:=`(
      vcf_id = NA_character_,
      Status = "Exclude"
    )
  ][
    i = is.na(vcf_id),
    j = is_related := NA
  ]

  merge(
    x = dt_snp,
    y =   methy_sample_sheet[
      j = .SD,
      .SDcols = grep("^CellT_|^Sample_|^Sentrix_|qc_sex_discrepancy|call_rate", names(methy_sample_sheet))
    ], # Add cell components and QC metrics
    by = "Sample_Name",
    all.x = TRUE
  )[
    i = (qc_sex_discrepancy | is.na(call_rate)),
    j = Status := "Exclude"
  ][
    j = if (.N > 1) .SD[!Status %in% "Exclude" & call_rate == max(call_rate, na.rm = TRUE)] else .SD,
    by = "Sample_Name"
  ]
}

#' do_meqtl
#' @import data.table
#' @import IlluminaHumanMethylationEPICanno.ilm10b5.hg38
#' @import future_apply
#' @import utils
do_meqtl <- function(
  phenotype,
  model,
  beta_file,
  vcfs,
  vep,
  path,
  epic_annot_pkg = "IlluminaHumanMethylationEPICanno.ilm10b5.hg38",
  bin_path = list(
    qtltools = "QTLtools",
    bcftools =  "/usr/bin/bcftools",
    tabix = "/usr/bin/tabix",
    bgzip = "/usr/bin/bgzip"
  ),
  cis_window = 500000,
  n_chunk = 20
) {
  tmp_dirs <- (function(x) `names<-`(file.path(path, x), x))(
    c("methylation", "genotypes", "covariates", "qtl", "qtl_annotated", "qtl_combined")
  )
  invisible(sapply(tmp_dirs, dir.create, recursive = TRUE, showWarnings = FALSE, mode = "0775"))

  beta_matrix <- (function(x) log2(x) - log2(1 - x))(
    as.matrix(data.table::fread(file = beta_file, header = TRUE), "cpg_id")
  )

  epic_qc_annot <- get(utils::data("Locations", package = epic_annot_pkg))
  epic_qc_annot <- epic_qc_annot[intersect(rownames(beta_matrix), rownames(epic_qc_annot)), ]
  epic_qc_annot <- na.exclude(as.data.frame(epic_qc_annot))
  epic_qc_annot <- epic_qc_annot[epic_qc_annot[["chr"]] %in% sprintf("chr%d", 1:22), ]
  epic_qc_annot[["chr"]] <- as.numeric(gsub("chr", "", epic_qc_annot[["chr"]]))

  epic_qc <- beta_matrix[
    rownames(epic_qc_annot),
    as.character(phenotype[["Sample_ID"]])
  ]
  colnames(epic_qc) <- phenotype[["vcf_id"]]

  devnull <- merge(
    x = data.table::as.data.table(epic_qc, keep.rownames = "CpG"),
    y = data.table::as.data.table(epic_qc_annot, keep.rownames = "CpG"),
    by = "CpG"
  )[
    j = (c("#Chr", "start", "end", "grp")) := list(chr, pos, pos, chr)
  ][
    j = .SD,
    .SDcols = c("grp", "#Chr", "start", "end", "CpG", phenotype[["vcf_id"]])
  ][
    order(`#Chr`, start)
  ][
    j = (function(x, y) {
      file <- sprintf("%s/chr%02d.bed", tmp_dirs[["methylation"]], unique(y))
      data.table::fwrite(x, file = file, quote = FALSE, sep = "\t")
      system(sprintf("%s -f %s", bin_path[["bgzip"]], file))
      system(sprintf("%s -p bed -f %s.gz", bin_path[["tabix"]], file))
      return(TRUE)
    })(.SD, `#Chr`),
    by = "grp"
  ]

  file_con <- gzfile(file.path(tempdir(), "nominal_header.txt.gz"), "w")
  cat(
    c("cpg_id", "rs_id", "distance_bp", "pvalue", "slope\n"),
    sep = " ",
    file = file_con
  )
  close(file_con)

  file_con <- gzfile(file.path(tempdir(), "permutation_header.txt.gz"))
  cat(
    c(
      "cpg_id", "variants_cis", "mle_shape1_beta", "mle_shape2_beta",
      "dummy", "best_rs_id", "distance_bp", "pvalue", "slope",
      "permutation_pvalue", "downstream_pvalue\n"
    ),
    sep = " ",
    file = file_con
  )
  close(file_con)

  data.table::fwrite(
    x = phenotype[j = .SD, .SDcols = "vcf_id"],
    file = sprintf("%s/keep.samples", tmp_dirs[["genotypes"]]),
    sep = " ",
    col.names = FALSE
  )

  list_vcfs <- future.apply::future_sapply(
    X = vcfs[sub("\\.vcf.gz", "", basename(vcfs)) %in% 1:22],
    vep_file = vep,
    bin_path = bin_path,
    output_genotypes = tmpdirs[["genotypes"]],
    future.globals = FALSE,
    future.packages = "data.table",
    FUN = function(vcf, vep_file, output_genotypes, bin_path) {
      vcf_file <- sprintf("%s/%s", output_genotypes, basename(vcf))
      system(paste(
        bin_path[["bcftools"]],
          "+fill-tags", vcf,
       "|",
        bin_path[["bcftools"]],
          "view",
          "--min-af 0.05",
          "--exclude 'INFO/INFO < 0.8'",
          "--min-alleles 2 --max-alleles 2 --types snps",
          "--force-samples",
          # "--no-update",
          "--samples-file", sprintf("%s/keep.samples", output_genotypes),
       "|",
        bin_path[["bcftools"]],
          "annotate",
          "--annotations", vep_file,
          "--header-lines", sub("_formatted.tsv.gz", ".header", vep_file),
          "--columns CHROM,POS,Gene,Symbol,rsid",
          "--set-id '%INFO/rsid'",
        "|",
        bin_path[["bcftools"]],
          "annotate",
          "--set-id +'%CHROM:%POS:%REF:%ALT'",
          "--output-type z --output", vcf_file
      ))
      vcf_file
    }
  )
  names(list_vcfs) <- sub("\\.vcf.gz", "", basename(list_vcfs))

  for (ichr in sprintf("chr%02d", 1:22)) {
    local({
      for (ianalysis in c("nominal", "permutation")) {
        future_apply::future_lapply(
          X = seq_len(n_chunk),
          future.globals = FALSE,
          ivcf = list_vcfs[ichr],
          vep_file = vep,
          tmp_dirs = tmp_dirs,
          cis_window = cis_window,
          ianalysis = ianalysis,
          ichr = ichr,
          n_chunk = n_chunk,
          cis_window = cis_window,
          FUN = function(ichunk, ivcf, tmp_dirs, cis_window, ianalysis, ichr, n_chunk, bin_path) {
            system(paste(bin_path[["qtltools"]],
              "cis",
              "--silent",
              "--vcf", ivcf,
              "--bed", sprintf("%s/%s.bed.gz", tmp_dirs[["methylation"]], ichr),
              "--cov", sprintf("%s/covariates.txt.gz", tmp_dirs[["covariates"]]),
              "--window", cis_window,
              ifelse(ianalysis == "nominal", "--nominal 1", "--permute 1000"),
              "--chunk", ichunk, n_chunk,
              "--out", sprintf("%s/%s_%s_%03d.txt.gz", tmp_dirs[["qtl"]], ianalysis, ichr, ichunk)
            ))
          }
        )

        system(sprintf(
          "zcat %s/%s_header.txt.gz %s/%s_%s_*.txt.gz | %s -c > %s/%s_%s.txt.gz",
          tempdir(), ianalysis,
          tmp_dirs[["qtl"]], ianalysis, ichr,
          bin_path[["bgzip"]],
          tmp_dirs[["qtl"]], ianalysis, ichr
        ))

        unlink(list.files(
          path = output_qtl,
          pattern = sprintf("%s_%s_[0-9]*.txt.gz", ianalysis, ichr),
          full.names = TRUE
        ))

        data.table::fwrite(
          x = Reduce(
            f = function(x, y) merge(x, y, by = "cpg_id", all.x = TRUE),
            x = list(
              data.table::fread(sprintf("%s/%s_%s_%s.txt.gz", tmp_dirs[["qtl"]], project_name, ianalysis, ichr)),
              data.table::as.data.table(epic_qc_annot, keep.rownames = "cpg_id"),
              data.table::as.data.table(get(utils::data("Other", package = epic_annot_pkg)), keep.rownames = "cpg_id")[
                j = list(
                  UCSC_RefGene_Name = paste(unique(tstrsplit(data.table::UCSC_RefGene_Name, split = ";")), collapse = ";")
                ),
                by = "cpg_id"
              ]
            )
          ),
          file = sprintf("%s/%s_%s_%s.txt.gz", tmp_dirs[["qtl_annotated"]], project_name, ianalysis, ichr)
        )
      }
    })
  }

  for (ianalysis in c("nominal", "permutation")) {
    local({
      data.table::fwrite(
        x = data.table::rbindlist(future_lapply(
          X = sprintf("chr%02d", 1:22),
          future.globals = FALSE,
          future.packages = "data.table",
          tmp_dirs = tmp_dirs[["qtl_annotated"]],
          project_name = project_name,
          ianalysis = ianalysis,
          FUN = function(ichr, tmp_dirs, project_name, ianalysis) {
            data.table::fread(sprintf("%s/%s_%s_%s.txt.gz", tmp_dirs[["qtl_annotated"]], project_name, ianalysis, ichr))
          }
        )),
        file = sprintf("%s/%s_%s.txt.gz", tmp_dirs[["qtl_combined"]], project_name, ianalysis)
      )
    })
  }

  on.exit(
    unlink(
      x = c(
        gzfile(file.path(tempdir(), "nominal_header.txt.gz")),
        gzfile(file.path(tempdir(), "permutation_header.txt.gz")),
        setdiff(
          list.files(path, full.names = TRUE),
          list.files(path, pattern = "qtl_combined", full.names = TRUE)
        )
      ),
      recursive = TRUE,
      force = TRUE
    )
  )

  list.files(path, pattern = "qtl_combined", full.names = TRUE)
}
