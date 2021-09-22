### global libraries ===============================================================================
library(targets)
library(tarchetypes)
library(here)
library(data.table)
library(future)
library(future.callr)

# tar_option_set(cue = tar_cue(mode = "never"))

# targets::tar_renv(extras = "visNetwork", path = "scripts/_dependencies.R")


### project setup ==================================================================================
invisible(sapply(list.files(here("scripts"), pattern = "^tar-.*R$", full.names = TRUE), source, echo = FALSE))

plan(future.callr::callr, workers = 3)
# plan(multicore, workers = 40)
message(sprintf("Number of workers: %d", future::nbrOfWorkers()))
# setDTthreads(threads = 1)

tar_setup <- {list( # Setup project
  tar_target(project, gsub("(.*)_.*", "\\1", list.files(here(), pattern = ".Rproj$")), packages = "here"),
  tar_target(author, "Mickaël CANOUIL, *Ph.D.*"),
  tar_target(output_directory, here::here("outputs"), packages = "here")
)}


### targets ========================================================================================
tar_sample_sheet_qc <- {list(
  tar_target(ewas_sample_sheet_qc,
    command = qc_sample_sheet_ewas(
      phenotype = phenotypes,
      methylation = file.path(ma_export_directory, "EPIC_QC_phenotypes.csv")
    ),
    packages = "data.table"
  )
)}

tar_ewas <- {list(
  tar_target(ewas_models,
    command = tar_group(dplyr::group_by(
      data.frame(
        pretty_trait = c("Case/Control"),
        raw_trait = c("group"),
        covariates = c(
          paste(c("sex", "age", "bmi"), collapse = " + "),
          paste(c("sex", "age", "bmi", "cell"), collapse = " + ")
        )
      ),
      pretty_trait, raw_trait, covariates
    )),
    packages = "dplyr",
    iteration = "group"
  ),
  tar_target(ewas_results_file,
    command = do_ewas(
      data = ewas_sample_sheet_qc[!Status %in% "Exclude"], # phenotypes
      model = ewas_models,
      beta_file = file.path(ma_export_directory, "EPIC_QC_betavalues.csv.gz"),
      path = file.path(output_directory, "ewas"),
      epic_annot_pkg = "IlluminaHumanMethylationEPICanno.ilm10b5.hg38"
    ),
    pattern = map(ewas_models),
    iteration = "list",
    packages = c(
      "here", "data.table", "stats", "utils", "future.apply", "limma",
      "IlluminaHumanMethylationEPICanno.ilm10b5.hg38"
    ),
    format = "file"
  ),
  tar_target(ewas_results_manhattan,
    command = plot_manhattan_ewas(file = ewas_results_file, model = ewas_models),
    pattern = map(ewas_models, ewas_results_file),
    iteration = "list",
    packages = c("ggplot2", "ggtext", "data.table", "ggrepel", "scales")
  ),
  tar_target(ewas_results_volcano,
    command = plot_volcano_ewas(file = ewas_results_file, model = ewas_models),
    pattern = map(ewas_models, ewas_results_file),
    iteration = "list",
    packages = c("ggplot2", "ggtext", "data.table", "ggrepel", "scales")
  ),
  tar_target(ewas_results_pp,
    command = plot_pp_ewas(file = ewas_results_file, model = ewas_models),
    pattern = map(ewas_models, ewas_results_file),
    iteration = "list",
    packages = c("ggplot2", "ggtext", "data.table", "stats")
  ),
  tar_render(ewas_report,
    path = here("scripts/slides/ewas_report.Rmd"),
    output_dir = here("reports"),
    packages = c(
      "xaringan",
      "here", "knitr", "ragg", "ggplot2", "ggtext",
      "patchwork", "data.table", "gt", "scales",
      "showtext", "svglite", "katex",
      "targets", "bacon", "utils"
    )
  )
)}

list(
  tar_setup,
  tar_sample_sheet_qc,
  tar_ewas
)
