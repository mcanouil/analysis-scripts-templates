### global libraries ===============================================================================
library(targets)
library(tarchetypes)
library(here)
library(data.table)
# library(future)
# library(future.callr)

# targets::tar_renv(extras = "visNetwork", path = "scripts/_dependencies.R")
# renv::install("gabraham/flashpca/flashpcaR")
# renv::install(c("bioc::S4Vectors", "bioc::biomaRt", "bioc::tximport", "bioc::DESeq2", "bioc::MatrixGenerics"))
# renv::hydrate()


### project setup ==================================================================================
# Functions/scripts required: tar-twas.R, tar-biomart.R, tar-pval_trans.R
invisible(sapply(
  X = list.files(here("scripts", "tar-utils"), pattern = "^tar-.*R$", full.names = TRUE),
  FUN = source, echo = FALSE
))

# plan(future.callr::callr, workers = 40)
# plan(multicore, workers = 40)
# message(sprintf("Number of workers: %d", future::nbrOfWorkers()))
# setDTthreads(threads = 1)


### targets setup ==================================================================================
tar_setup <- list(
  tar_target(project, sub("(.*)_[^_]*\\.Rproj$", "\\1", list.files(here(), pattern = ".Rproj$")), packages = "here"),
  tar_target(author, "Mickaël Canouil, *Ph.D*"),
  tar_target(output_directory, here::here("outputs"), packages = "here")
)


### targets ========================================================================================
tar_sample_sheet_qc <- list(
  tar_target(twas_sample_sheet_qc,
    command = qc_sample_sheet_twas(
      phenotype = NULL,
      run_path = "/disks/RUN/Run_42/"
    ),
    packages = "data.table"
  )
)
tar_twas <- list(
  tar_target(twas_rna_level, c("ensembl_gene_id", "ensembl_transcript_id")[1]),
  tar_target(twas_models,
    command = tar_group(dplyr::group_by(
      data.frame(
        pretty_trait = c("Case/Control"),
        raw_trait = c("group"),
        covariates = c("")
      ),
      pretty_trait, raw_trait, covariates
    )),
    packages = "dplyr",
    iteration = "group"
  ),
  tar_target(twas_tximport,
    command = read_rsem(
      sample_sheet = twas_sample_sheet_qc,
      rna_level = twas_rna_level
    ),
    pattern = map(twas_rna_level),
    iteration = "list",
    packages = c("tximport", "readr")
  ),
  tar_target(twas_pca_plots,
    command = plot_pca_twas(
      txi = twas_tximport,
      sample_sheet = twas_sample_sheet_qc,
      pca_vars = c("group"),
      n_comp = 10,
      fig_n_comp = 3
    ),
    pattern = map(twas_tximport, twas_rna_level),
    iteration = "list",
    packages = c(
      "flashpcaR", "data.table", "ggplot2", "ggtext", "patchwork",
      "scales", "stats", "utils",
      "DESeq2", "MatrixGenerics"
    )
  ),
  tar_target(twas_biomart,
    command = get_biomart_information(
      ensembl_id = rownames(twas_tximport[["counts"]]),
      rna_level = twas_rna_level,
      organism = "hsapiens_gene_ensembl",
      version = unique(sub(
        pattern = ".*Ensembl-([0-9]+)\\.genes\\.results$",
        replacement = "\\1",
        x = basename(twas_sample_sheet_qc[["rnaseq_path"]])
      ))
    ),
    pattern = map(twas_rna_level, twas_tximport),
    iteration = "list",
    packages = c("biomaRt", "httr", "data.table")
  ),
  tar_target(twas_results_file,
    command = do_twas(
      txi = twas_tximport,
      sample_sheet = twas_sample_sheet_qc,
      model = twas_models,
      path = file.path(output_directory, "twas"),
      rna_level = twas_rna_level,
      biomart = twas_biomart
    ),
    pattern = cross(twas_models, map(twas_rna_level, twas_tximport, twas_biomart)),
    iteration = "list",
    packages = c("data.table", "DESeq2", "S4Vectors", "MatrixGenerics", "utils", "stats"),
    format = "file"
  ),
  tar_target(twas_results_volcano,
    command = plot_volcano_twas(file = twas_results_file, model = twas_models),
    pattern = map(twas_models, twas_results_file),
    iteration = "list",
    packages = c("ggplot2", "ggtext", "data.table", "ggrepel", "scales")
  ),
  tar_render(twas_report,
    path = here("slides/twas_report.Rmd"),
    output_dir = here("reports"),
    packages = c(
      "xaringan",
      "here", "knitr", "ragg", "ggplot2", "ggtext",
      "patchwork", "data.table", "gt", "scales",
      "showtext", "svglite", "katex",
      "targets", "bacon", "utils"
    )
  )
)

list(
  tar_setup,
  tar_sample_sheet_qc,
  tar_twas
)
