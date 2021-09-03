#' do_gwas
#' @import data.table
#' @import stats
#' @import future.apply future_lapply
do_gwas <- function(
  data,
  model,
  vcfs,
  vep,
  path,
  bin_path = list(
    bcftools = "/usr/bin/bcftools",
    plink2 = "plink2"
  )
) {
  tmpdir <- file.path(tempdir(), "gwas_plink2", model[["raw_trait"]])
  dir.create(path = tmpdir, recursive = TRUE, mode = "0777")
  on.exit(unlink(tmpdir, recursive = TRUE))

  vcfs <- vcfs[sub("\\.vcf.gz", "", basename(vcfs)) %in% 1:22]

  covariates <- all.vars(as.formula(sprintf(" ~ %s", model[["covariates"]])))
  dt <- unique(data.table::setnames(
    x = data.table::as.data.table(data)[
      j = na.exclude(.SD),
      .SDcols = c("vcf_id", model[["raw_trait"]], covariates)
    ],
    old = "vcf_id",
    new = "#IID"
  ))

  if (length(sex_covariate <- grep("^sex", covariates, value = TRUE)) > 1) {
    stop(sprintf(
      'Only one column containing "sex" can exist in the model, not %s: "%s"',
      length(sex_covariate),
      paste(sex_covariate, collapse = '", "')
    ))
  }

  model_type <- if (data.table::uniqueN(dt[[model[["raw_trait"]]]]) == 2) {
    "logistic"
  } else {
    "linear"
  }
  basename_file <- sprintf("%s/%s_%s", tmpdir, model_type, model[["raw_trait"]])

  data.table::fwrite(
    x = dt[j = .SD, .SDcols = "#IID"],
    file = sprintf("%s.samples", basename_file),
    sep = " "
  )

  data.table::fwrite(
    x = dt[j = .SD, .SDcols = c("#IID", model[["raw_trait"]])],
    file = sprintf("%s.pheno", basename_file),
    sep = " "
  )

  data.table::fwrite(
    x = dt[j = .SD, .SDcols = c("#IID", setdiff(covariates, sex_covariate))],
    file = sprintf("%s.cov", basename_file),
    sep = " "
  )

  data.table::fwrite(
    x = data.table::setnames(
      x = data.table::copy(dt)[j = .SD, .SDcols = c("#IID", sex_covariate)],
      old = sex_covariate,
      new = "SEX"
    ),
    file = sprintf("%s.sex", basename_file),
    sep = " "
  )

  message("Formatting VCFs and performing PLINK2 regression ...")

  list_results <- future.apply::future_lapply(
    X = vcfs,
    basename_file = basename_file,
    vep_file = vep,
    bin_path = bin_path,
    future.globals = FALSE,
    future.packages = "data.table",
    FUN = function(vcf, basename_file, vep_file, bin_path) {
      vcf_file <- sprintf("%s__%s", basename_file, basename(vcf))
      results_file <- sub("\\.vcf.gz", "", vcf_file)

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
          "--samples-file", sprintf("%s.samples", basename_file),
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

      system(paste(
        bin_path[["plink2"]],
        "--vcf", vcf_file, "dosage=DS",
        "--mach-r2-filter",
        "--threads 120",
        "--glm sex",
        "--keep", sprintf("%s.samples", basename_file),
        "--update-sex", sprintf("%s.sex", basename_file),
        "--pheno", sprintf("%s.pheno", basename_file),
        "--covar", sprintf("%s.cov", basename_file),
        "--covar-variance-standardize",
        "--silent",
        "--out", results_file
      ))

      annot <- data.table::fread(
        cmd = paste(bin_path[["bcftools"]], "view --drop-genotypes", vcf_file),
        skip = "#CHROM"
      )
      annot <- annot[
        j = list(
          .SD,
          data.table::rbindlist(
            l = lapply(
              X = strsplit(INFO, ";"),
              FUN = function(x) {
                all_fields <- strsplit(x, "=")
                out <- data.table::transpose(all_fields[sapply(all_fields, length) > 1])
                data.table::setnames(
                  x = data.table::as.data.table(do.call("rbind", out[-1])),
                  old = out[[1]]
                )
              }
            ),
            use.names = TRUE,
            fill = TRUE
          )[
            j = lapply(.SD, function(x) {
              xout <- as.character(x)
              data.table::fifelse(
                test = xout %in% c(".", "-"),
                yes = NA_character_,
                no = xout
              )
            })
          ]
        ),
        .SDcols = !c("INFO", "QUAL", "FILTER")
      ]
      data.table::setnames(
        x = annot,
        old = function(x) sub("^\\.SD\\.\\.*", "", x)
      )

      results <- data.table::setnames(
        x = data.table::fread(
          file = list.files(
            path = dirname(results_file),
            pattern = sprintf("%s\\..*\\.glm\\..*", basename(results_file)),
            full.names = TRUE
          )
        ),
        old = function(x) sub("^#", "", x)
      )

      data.table::fwrite(
        x = data.table::merge.data.table(
          x = results[TEST %in% "ADD" & !is.na(P), -c("TEST")],
          y = annot,
          by = c("CHROM", "POS", "ID", "REF", "ALT"), # intersect(names(results), names(annot))
        )[MAF >= 0.05, .SDcols = !c("QUAL", "FILTER")],
        file = sprintf("%s.results.gz", results_file)
      )

      sprintf("%s.results.gz", results_file)
    }
  )

  results_file <- sprintf("%s/gwas_%s.csv.gz", path, model[["raw_trait"]])
  dir.create(
    path = dirname(results_file),
    recursive = TRUE,
    mode = "0755",
    showWarnings = FALSE
  )

  message("Aggregating PLINK2 results ...")

  data.table::fwrite(
    x = data.table::setcolorder(
      x = data.table::rbindlist(lapply(list_results, data.table::fread), use.names = TRUE)[
        j = `:=`(
          FDR = stats::p.adjust(P, method = "BH"),
          Bonferroni = stats::p.adjust(P, method = "bonferroni"),
          trait = model[["pretty_trait"]],
          covariates = model[["covariates"]]
        )
      ][order(P)],
      neworder = c("trait", "covariates")
    ),
    file = results_file
  )

  message(sprintf('Writing results to "%s"!', results_file))

  results_file
}

#' plot_manhattan_gwas
#' @import ggplot2
#' @import ggtext
#' @import data.table
#' @import ggrepel
plot_manhattan_gwas <- function(file, model) {
  dt <- data.table::fread(file)[
    i = P <= 0.05 / .N,
    j = gene_label := data.table::fifelse(Symbol == "", NA_character_, Symbol)
  ][
    j = gene_label_min := data.table::fifelse(
      test = P == min(P, na.rm = TRUE),
      yes = gsub(",", ";", gene_label),
      no = NA_character_
    ),
    by = "gene_label"
  ][
    i = P > 0.05,
    j = P := NA_real_
  ][
    j = file := basename(file)
  ][
    j = c("CHROM", "POS", "P", "gene_label_min", "gene_label")
  ]

  draw_manhattan(
    data = dt,
    x = "POS",
    y = "P",
    chr = "CHROM",
    label_y = "P-value",
    alpha = 0.05 / nrow(dt)
  ) +
    ggrepel::geom_label_repel(
      mapping = ggplot2::aes(label = gene_label_min),
      stat = "manhattan",
      show.legend = FALSE,
      min.segment.length = 0,
      # direction = "x",
      size = 1.75
    ) +
    ggplot2::labs(
      title = model[["pretty_trait"]],
      subtitle = toupper(paste(model[["group"]], "= ", model[["covariates"]]))
    ) +
    ggplot2::theme_minimal(base_family = "Verdana") +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggtext::element_markdown(),
      plot.subtitle = ggtext::element_markdown(face = "italic"),
      axis.text.y = ggtext::element_markdown(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.position = "none"
    )
}
