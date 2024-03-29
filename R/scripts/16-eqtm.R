message(timestamp(quiet = TRUE))
### Project Setup ==================================================================================
library(here)
project_name <- sub("(.*)_*\\..*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "16-eqtm")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

debug <- FALSE # FALSE or number of chunks to run

workers <- 40
workers_multiplier <- 5
set.seed(20210727)

cis_window <- 1e6 # in base pair

species <- "hsapiens_gene_ensembl"
build <- 38
version <- 104

data_directory <- file.path("/disks/DATA/Projects", project_name, "QC")
run_directory <- "/disks/RUN/Run_XXX/Output/RSEM"

phenotype <- here("docs", "phenotype.xlsx")

do_rna_level <- "isoforms" #  c("isoforms", "genes")
default_covariates <- "sex + age + bmi + PC01 + PC02"


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(data.table)

  library(future)
  library(future.apply)

  library(tximport)
  library(DESeq2)
  library(matrixStats)

  library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
  library(biomaRt)

  library(broom)
  library(qvalue)

  library(bench)
})

message("##------           Setup          ------##")

##------ future ------##
workers <- min(workers, availableCores())
if (interactive()) plan(sequential) else plan(multicore, workers = workers)
message(sprintf("Number of workers: %d", workers))

##------ data.table ------##
setDTthreads(threads = workers)


### Setup biomaRt ==================================================================================
get_mart <- quote(useEnsembl(
  biomart = "ensembl",
  dataset = species,
  version = version,
  GRCh = if (build == 37) build else NULL
))

mart <- try(eval(get_mart), silent = TRUE)
if (inherits(mart, "try-error")) mart <- eval(get_mart)
ensembl_version <- sprintf("GRCh%d-%d", build, version)


### Analysis =======================================================================================
message("##------          Format          ------##")
#### Phenotypes ------------------------------------------------------------------------------------
phenotype_matrix <- fread(
  file = file.path(data_directory, "EPIC", "EPIC_QC_phenotypes.csv"),
  colClasses = c("Sample_ID" = "character")
)

sample_sheet_qc <- merge(
  x = fread(
    file = phenotype,
    colClasses = c("Sample_ID" = "character")
  ),
  y = phenotype_matrix[
    j = .SD,
    .SDcols = grep("CellT|Sample_ID|ID_unique_interv", colnames(phenotype_matrix), value = TRUE)
  ],
  by.x = "Sample_ID",
  by.y = "Sample_ID"
)[
  j = .SD,
  .SDcols = c("Sample_ID", unlist(strsplit(default_covariates, " \\+ ")))
]

#### EPIC data -------------------------------------------------------------------------------------
raw_beta_matrix <- fread(
  file = file.path(data_directory, "EPIC", "EPIC_QC_betavalues.csv.gz"),
  header = TRUE,
  select = c("cpg_id", sample_sheet_qc[["Sample_ID"]])
)

methyl_annot <- suppressWarnings(as.data.table( # In .local(x, row.names, optional, ...) : Arguments in '...' ignored
  x = IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other[,
    c(
      "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group",
      "CHR_hg38", "Start_hg38", "End_hg38", "Strand_hg38"
    )
  ],
  keep.rownames = "cpg_id"
))[
  i = !is.na(Start_hg38) & cpg_id %in% intersect(raw_beta_matrix[["cpg_id"]], cpg_id),
  j = list(CHR_hg38, Start_hg38 = as.integer(Start_hg38), cpg_id)
]

raw_beta_matrix <- (function(x) log2(x) - log2(1 - x))(
  as.matrix(raw_beta_matrix[cpg_id %in% methyl_annot[["cpg_id"]]], "cpg_id")
)

message(sprintf(
  "Number of CpGs: %s",
  format(nrow(raw_beta_matrix), big.mark = ",")
))

#### RNA-seq data & eQTM ---------------------------------------------------------------------------
for (rna_level in do_rna_level) {
  rna_level_name <- unname(c("genes" = "gene", "isoforms" = "transcript")[rna_level])

  txi_counts <- suppressMessages(tximport(
    files = setNames(
      object = sprintf(
        "%s/sample%04d_Ensembl-%s.%s.results",
        run_directory, as.numeric(sample_sheet_qc[["Sample_ID"]]),
        sub(".*-", "", ensembl_version), rna_level
      ),
      nm = sample_sheet_qc[["Sample_ID"]]
    ),
    type = "rsem",
    txIn = rna_level == "isoforms",
    txOut = rna_level == "isoforms",
    countsFromAbundance = "no",
  ))
  txi_counts$length[txi_counts$length == 0] <- 1

  dds <- suppressMessages(DESeqDataSetFromTximport(
    txi = txi_counts,
    colData = data.frame(intercept = rep(1, ncol(txi_counts[["counts"]]))),
    design = ~ 1
  ))
  counts_vst <- suppressMessages(assay(# filter on Raw counts and normalise with transcript/library length
    vst(dds[rowVars(counts(dds)) != 0 & rowMeans(counts(dds)) > 1 & rowMedians(counts(dds)) > 0, ])
  ))

  ensembl_ids <- setnames(
    data.table(
      id = rownames(counts_vst),
      id_noversion = sub("\\..*$", "", rownames(counts_vst))
    ),
    function(x) sprintf("%s_%s", rna_level_name, x)
  )

  annot_biomart <- setnames(
    x = setDT(getBM(
      attributes = c(
        sprintf("ensembl_%s_id", rna_level_name),
        "chromosome_name",
        list(
          gene = c("start_position", "end_position"),
          transcript = c("transcript_start", "transcript_end")
        )[[rna_level_name]],
        "external_gene_name"
      ),
      filters = sprintf("ensembl_%s_id", rna_level_name),
      values = list(unique(ensembl_ids[[sprintf("%s_id_noversion", rna_level_name)]])),
      mart = mart
    )),
    old = function(x) sub("^transcript_|_position$", "", x)
  )[
    i = order(chromosome_name, start),
    j = `:=`(
      chr = sprintf("chr%s", chromosome_name),
      start_cis_window = fifelse(
        test = start < cis_window,
        yes = 0,
        no  = start - cis_window
      ),
      end_cis_window = start + cis_window
    )
  ]

  rna_annot <- merge(
    x = ensembl_ids,
    y = annot_biomart,
    by.x = sprintf("%s_id_noversion", rna_level_name),
    by.y = sprintf("ensembl_%s_id", rna_level_name)
  )[j = ensembl := ensembl_version]

  message(sprintf(
    "Number of %ss: %s",
    rna_level_name,
    format(nrow(rna_annot), big.mark = ",")
  ))


  #### Make Pairs ----------------------------------------------------------------------------------
  rna_annot <- setkey(rna_annot, chr, start_cis_window, end_cis_window)[chr %in% sprintf("chr%d", 1:22)]
  methyl_annot <- setkey(methyl_annot, CHR_hg38, Start_hg38)[
    CHR_hg38 %in% sprintf("chr%d", 1:22)
  ][
    j = position_cpg := Start_hg38
  ]

  cis_cpg_gene_pairs_info <- methyl_annot[
    i = rna_annot,
    j = list(count_id = get(sprintf("%s_id", rna_level_name)), start, cpg_id, position_cpg),
    on = list(CHR_hg38 = chr, Start_hg38 >= start_cis_window, Start_hg38 <= end_cis_window),
    by = .EACHI,
    nomatch = NULL
  ][
    j = `:=`(
      dist = start - position_cpg,
      W = as.integer(factor(count_id)) %/% (workers * workers_multiplier)
    )
  ][
    i = order(count_id, cpg_id),
    j = c("count_id", "start", "cpg_id", "CHR_hg38", "position_cpg", "dist", "W")
  ]

  cis_cpg_gene_pairs <- cis_cpg_gene_pairs_info#[j = c("count_id", "cpg_id", "W")] # keep or not position/distance/etc.

  counts_vst <- counts_vst[unique(cis_cpg_gene_pairs[["count_id"]]), , drop = FALSE]
  beta_matrix <- raw_beta_matrix[unique(cis_cpg_gene_pairs[["cpg_id"]]), , drop = FALSE]

  #### eQTM ----------------------------------------------------------------------------------------
  message(timestamp(quiet = TRUE))
  message("##------           eQTM           ------##")

  message(sprintf(
    "Number of cis-pairs: %s",
    format(nrow(cis_cpg_gene_pairs), big.mark = ",")
  ))
  message(sprintf(
    "Number of %ss considered in the cis-pairs: %s",
    rna_level_name,
    format(nrow(counts_vst), big.mark = ",")
  ))
  message(sprintf(
    "Number of CpGs considered in the cis-pairs: %s",
    format(nrow(beta_matrix), big.mark = ",")
  ))
  message(sprintf(
    "Number of chunks of %ss: %s",
    rna_level_name,
    format(max(cis_cpg_gene_pairs[["W"]]), big.mark = ",")
  ))

  ##------ DEBUGGING ------##
  if (!isFALSE(debug)) cis_cpg_gene_pairs <- cis_cpg_gene_pairs[W %in% 0:(debug - 1)]
  ##------    END    ------##

  message(sprintf("Time taken for eQTM: %s", bench_time({
    results <- cis_cpg_gene_pairs[j = list(W = unique(W), file = NA_character_), by = "count_id"]

    progress_steps <- unique(ceiling(c(
      min(cis_cpg_gene_pairs[["W"]]),
      seq(
        from = min(cis_cpg_gene_pairs[["W"]]),
        to = max(cis_cpg_gene_pairs[["W"]]),
        length.out = 100
      ),
      max(cis_cpg_gene_pairs[["W"]])
    )))

    for (ichunk in unique(cis_cpg_gene_pairs[["W"]])) {
      setDTthreads(threads = workers)
      cis_cpg_gene_pairs_data <- cis_cpg_gene_pairs[W %in% ichunk][
        j = list(
          data = list(
            as.data.table(
              x = t(rbind(
                counts_vst[count_id, , drop = FALSE],
                beta_matrix[cpg_id, , drop = FALSE]
              )),
              keep.rownames = "ABOS_ID"
            )
          )
        ),
        by = "count_id"
      ]
      setDTthreads(threads = 1)

      chunk_data <- cis_cpg_gene_pairs_data[["data"]]

      files <- future_sapply(
        X = chunk_data,
        cov_formula = default_covariates,
        pheno = sample_sheet_qc,
        future.packages = c("data.table", "stats", "broom", "progress"),
        future.globals = FALSE,
        future.chunk.size = workers_multiplier,
        FUN = function(data, cov_formula, pheno) {
          setDTthreads(threads = 1)

          lm_eqtm <- function(formula, data) {
            setDT(tidy(lm(
              formula = as.formula(sprintf("%s ~ mvalue + %s", names(data)[[2]], formula)),
              data = data
            )))[
              term %in% "mvalue"
            ][
              j = `:=`(response = names(data)[[2]], term = NULL)
            ]
          }

          file <- sprintf("%s/%s.csv.gz", tempdir(), names(data)[[2]])

          fwrite(
            x = merge(
              x = melt(data = data, id.vars = 1:2, variable.name = "cpg_id", value.name = "mvalue"),
              y = pheno,
              by = "Sample_ID"
            )[
              j = lm_eqtm(formula = cov_formula, data = .SD),
              by = "cpg_id"
            ],
            file = file
          )

          file
        }
      )

      results[W %in% ichunk, file := files]

      if (ichunk %in% progress_steps) {
        message(sprintf("  - Completion = %0.2f %%", ichunk / max(progress_steps) * 100))
      }
    }
  })[[2]]))

  message(sprintf("Completed at: %.2f %%", results[j = sum(!is.na(file)) / .N * 100]))

  output_file <- sprintf("%s/%s_eQTM_%s.csv.gz", output_directory, project_name, rna_level_name)

  fwrite(
    x = merge(
      x = setnames(cis_cpg_gene_pairs[j = -c("W")], "count_id", sprintf("%s_id", rna_level_name)),
      y = setnames(rbindlist(lapply(results[["file"]], fread)), "response", sprintf("%s_id", rna_level_name)),
      by = c(sprintf("%s_id", rna_level_name), "cpg_id")
    )[j = fdr_storey := qvalue(p.value)$qvalues],
    file = output_file
  )
  unlink(x = results[["file"]])

  message(sprintf('Results available in: "%s"', sub(here(), ".", output_file)))
}


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
message(timestamp(quiet = TRUE))
