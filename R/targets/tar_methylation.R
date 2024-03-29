### global libraries ===============================================================================
library(targets)
library(tarchetypes)
library(here)
library(data.table)
# library(future)
# library(future.callr)

# targets::tar_renv(extras = c("visNetwork", "umr1283"), path = "scripts/_dependencies.R")

# options(configure.args = "--disable-threading"); renv::install("bmbolstad/preprocessCore", force = TRUE)
# renv::hydrate()
# renv::install(c(
#   "gabraham/flashpca/flashpcaR",
#   "RefFreeEWAS@2.2",
#   "bioc::sva",
#   "achilleasNP/IlluminaHumanMethylationEPICmanifest",
#   "bioc::ENmix",
#   "bioc::minfi",
#   "bioc::IlluminaHumanMethylation450kmanifest",
#   "bioc::ChAMPdata",
#   "achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38",
#   "bioc::FlowSorted.CordBloodCombined.450k",
#   "bioc::illuminaio",
#   "bioc::FlowSorted.Blood.EPIC"
# ))


### project setup ==================================================================================
# Functions/scripts required: tar-qc_idats.R, tar-pval_trans.R, tar-clean_sample_sheet.R
invisible(sapply(
  X = list.files(here("scripts", "tar-utils"), pattern = "^tar-.*R$", full.names = TRUE),
  FUN = source, echo = FALSE
))

# plan(future.callr::callr, workers = 40)
# plan(multicore, workers = 40)
# message(sprintf("Number of workers: %d", future::nbrOfWorkers()))
# setDTthreads(threads = 1)


### targets setup ==================================================================================
tar_setup <- list( # Setup project
  tar_target(project, sub("(.*)_[^_]*\\.Rproj$", "\\1", list.files(here(), pattern = ".Rproj$")), packages = "here"),
  tar_target(author, "Mickaël CANOUIL, *Ph.D.*"),
  tar_target(qc_data_path, "/disks/DATA/Projects"),
  tar_target(ma_run, "/disks/RUN/Array/Run/Results/epic/"),
  tar_target(ma_sample_sheet,
    command = file.path(ma_run, "samplesheet.csv"),
    format = "file"
  ),
  tar_target(ma_csv,
    command = clean_sample_sheet(sample_sheet = ma_sample_sheet, design = NULL),
    format = "file",
    packages = c("data.table", "readxl", "here")
  )
)


### targets ========================================================================================
tar_methylation <- list( # Methylation Array (ma)
  tar_target(ma_params,
    command = list(
      output_directory = here("outputs", "qc_ma"),
      csv_file = ma_csv,
      data_directory = ma_run,
      array = "Illumina EPIC",
      annotation = "ilm10b5.hg38",
      filter_snps = FALSE,
      filter_non_cpg = TRUE,
      filter_xy = TRUE,
      filter_multihit = TRUE,
      filter_beads = TRUE,
      population = NULL,
      bead_cutoff = 0.05,
      detection_pvalues = 0.01,
      filter_callrate = TRUE,
      callrate_samples = 0.95,
      callrate_probes = 1,
      sex_threshold = NULL,
      sex_colname = "sex_offspring",
      norm_background = "oob",
      norm_dye = "RELIC",
      norm_quantile = "quantile1",
      cell_tissue = "cordblood",
      pca_vars = c("Sample_Plate", "Sentrix_ID", "cohort"),
      max_labels = 15
    ),
    packages = "here"
  ),
  tar_target(ma_setup,
    command = dir.create(ma_params[["output_directory"]], showWarnings = FALSE, mode = "0775")
  ),
  tar_target(ma_data_idats,
    command = qc_idats(ma_params),
    packages = c(
      "data.table", "future", "future.apply",
      "IlluminaHumanMethylationEPICanno.ilm10b5.hg38",
      "IlluminaHumanMethylationEPICmanifest",
      # "IlluminaHumanMethylation450kmanifest",
      "ChAMPdata", "minfi", "utils", "ENmix", "illuminaio", "methods"
    )
  ),
  tar_target(ma_phenotypes,
    command = compute_phenotypes(
      data_mset = ma_data_idats[["mset"]],
      cell = ma_cell,
      sex_predicted = ma_check_sex
    ),
    packages = "data.table"
  ),
  tar_target(ma_cell,
    command = {
      estimate_cell_composition(
        data_rgset = ma_data_idats[["rgset"]],
        data_mset = ma_data_idats[["mset"]],
        cell_tissue = ma_params[["cell_tissue"]],
        array = sub(".* ", "", ma_params[["array"]]),
        n_cores = min(c(10, detectCores()))
      )
    },
    packages = c(
      "stats", "future", "future.apply",
      "minfi", "RefFreeEWAS",
      "FlowSorted.Blood.EPIC",
      "FlowSorted.CordBloodCombined.450k",
      "IlluminaHumanMethylation450kmanifest"
    )
  ),
  tar_target(ma_check_sex,
    command = check_sex(data_rgset = ma_data_idats[["rgset"]], sex_threshold = ma_sex_threshold),
    packages = c("minfi", "stats", "data.table")
  ),
  tar_target(ma_sex_threshold,
    command = compute_sex_threshold(
      data_rgset = ma_data_idats[["rgset"]],
      sex_threshold = ma_params[["sex_threshold"]]
    ),
    packages = c("minfi", "stats", "data.table")
  ),
  tar_target(ma_normalised_mset,
    command = normalise_mset(data_mset = ma_data_idats[["mset"]], phenotypes = ma_phenotypes),
    packages = c("minfi", "ENmix", "data.table", "sva")
  ),
  tar_target(ma_pca_mset_plots,
    command = mset_pca_plot(
      data = ma_data_idats,
      normalised_mset = ma_normalised_mset,
      pca_vars = ma_params[["pca_vars"]]
    ),
    packages = c(
      "flashpcaR", "data.table", "ggplot2", "ggtext", "patchwork",
      "scales", "stats", "utils"
    )
  ),
  tar_target(ma_callrate_plot,
    command = plot_callrate_ma(
      data = ma_phenotypes,
      callrate = ma_params[["callrate_samples"]],
      max_labels = ma_params[["max_labels"]]
    ),
    packages = c("data.table", "ggplot2", "ggrepel", "scales")
  ),
  tar_target(ma_sex_plot,
    command = plot_check_methylation_sex(
      data = ma_phenotypes,
      sex_threshold = ma_sex_threshold
    ),
    packages = c("data.table", "ggplot2", "ggrepel", "scales", "patchwork")
  ),
  tar_target(ma_cell_plot,
    command = plot_cell_composition(
      data = ma_phenotypes,
      max_labels = ma_params[["max_labels"]]
    ),
    packages = c("data.table", "ggplot2", "ggrepel", "scales", "patchwork", "ggdendro")
  ),
  tar_target(ma_export,
    command = export_ma_data(
      data_idats = ma_data_idats,
      mset = ma_normalised_mset,
      array = sub(".* ", "", ma_params[["array"]]),
      output_directory = ma_params[["output_directory"]]
    ),
    packages = c("readr", "data.table")
  ),
  tar_target(ma_export_directory,
    command = create_ma_export_directory(
      path = qc_data_path,
      project = project,
      array = sub(".* ", "", ma_params[["array"]])
    )
  ),
  tar_render(ma_qc_report,
    path = here("slides/methylation_array_qc_report.Rmd"),
    output_dir = here("reports"),
    params = ma_params,
    packages = c(
      "xaringan",
      "here", "knitr", "ragg", "ggplot2", "ggtext",
      "patchwork", "data.table", "gt", "scales",
      "targets"
    )
  ),
  tar_target(ma_save_qc,
    command = save_ma_qc_data(
      from = ma_export,
      report = here(grep("\\.html$", ma_qc_report, value = TRUE)),
      to = ma_export_directory
    ),
    packages = "here"
  )
)

list(
  tar_setup,
  tar_methylation
)
