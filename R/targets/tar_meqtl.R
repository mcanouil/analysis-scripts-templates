### global libraries ===============================================================================
library(targets)
library(tarchetypes)
library(here)
library(data.table)
library(future)
library(future.callr)

# targets::tar_renv(extras = "visNetwork", path = "scripts/_dependencies.R")


### project setup ==================================================================================
# Functions/scripts required: tar-meqtl.R,
invisible(sapply(list.files(here("scripts"), pattern = "^tar-.*R$", full.names = TRUE), source, echo = FALSE))

plan(future.callr::callr, workers = 3)
# plan(multicore, workers = 40)
message(sprintf("Number of workers: %d", future::nbrOfWorkers()))
# setDTthreads(threads = 1)


### targets setup ==================================================================================
tar_setup <- list( # Setup project
  tar_target(project, gsub("(.*)_.*", "\\1", list.files(here(), pattern = ".Rproj$")), packages = "here"),
  tar_target(author, "Mickaël CANOUIL, *Ph.D.*"),
  tar_target(output_directory, here::here("outputs"), packages = "here"),
  tar_target(bcftools, "/usr/bin/bcftools", format = "file"),
  tar_target(tabix, "/usr/bin/tabix", format = "file"),
  tar_target(bgzip, "/usr/bin/bgzip", format = "file"),
  tar_target(qtltools, "qtltools", format = "file")
)


### targets ========================================================================================
tar_sample_sheet_qc <- list(
  tar_target(meqtl_sample_sheet_qc,
    command = qc_sample_sheet_meqtl(
      phenotype = phenotypes,
      methylation = file.path(ma_export_directory, "EPIC_QC_phenotypes.csv"),
      exclusion = file.path(ga_export_directory, "quality-control-exclusion-checks.csv"),
      relatedness = file.path(ga_export_directory, "quality-control-relatedness.csv"),
      ethnicity = file.path(ga_export_directory, "quality-control-ethnicity.csv")
    ),
    packages = "data.table"
  )
)

tar_meqtl <- list(
  tar_target(meqtl_models,
    command = tar_group(dplyr::group_by(
      data.frame(
        pretty_trait = c("meQTL"),
        raw_trait = c("meQTL"),
        covariates = c(
          paste(c("sex", "age", "bmi"), collapse = " + "),
          paste(c("sex", "age", "bmi", "cell"), collapse = " + "),
          paste(c("sex", "age", "bmi", "cell", sprintf("PC%02d", 1:5)), collapse = " + ")
        )
      ),
      pretty_trait, raw_trait, covariates
    )),
    packages = "dplyr",
    iteration = "group"
  ),
  tar_target(meqtl_results_file,
    command = do_meqtl(
      phenotype = meqtl_sample_sheet_qc[!Status %in% "Exclude"], # phenotypes
      model = meqtl_models,
      beta_file = file.path(ma_export_directory, "EPIC_QC_betavalues.csv.gz"),
      vcfs = ga_imputed_vcf_crossmap,
      vep = veb_symbol_file,
      path = file.path(output_directory, "meqtl"),
      epic_annot_pkg = "IlluminaHumanMethylationEPICanno.ilm10b5.hg38",
      bin_path = list(
        qtltools = qtltools,
        bcftools =  bcftools,
        tabix = tabix,
        bgzip = bgzip
      ),
      cis_window = 500000,
      n_chunk = 20
    ),
    pattern = map(meqtl_models),
    iteration = "list",
    packages = c(
      "here", "data.table", "future.apply",
      "IlluminaHumanMethylationEPICanno.ilm10b5.hg38"
    ),
    format = "file"
  )
)

list(
  tar_setup,
  tar_sample_sheet_qc,
  tar_meqtl
)
