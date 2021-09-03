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
  tar_target(output_directory, here::here("outputs"), packages = "here"),
  tar_target(bcftools, "/usr/bin/bcftools", format = "file"),
  tar_target(plink2,
    command = download_plink2(
      url = "http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip",
      path = here("outputs", "plink2")
    ),
    format = "file",
    packages = "utils"
  )
)}


### targets ========================================================================================
tar_sample_sheet_qc <- list({
  tar_target(sample_sheet_qc,
    command = qc_sample_sheet(
      phenotype = harmonised_phenotypes,
      exclusion = file.path(ga_export_directory, "quality-control-exclusion-checks.csv"),
      relatedness = file.path(ga_export_directory, "quality-control-relatedness.csv"),
      ethnicity = file.path(ga_export_directory, "quality-control-ethnicity.csv"),
      methylation = ma_csv
    ),
    packages = "data.table"
  )
})

tar_gwas <- {list(
  tar_target(gwas_models,
    command = tar_group(dplyr::group_by(
      data.table(
        pretty_trait = "Control/Case",
        raw_trait = "group",
        covariates = paste(c(
          "sex", "age", "bmi",
          sprintf("PC%02d", 1:5)
        ), collapse = " + ")
      ),
      pretty_trait, raw_trait, covariates
    )),
    packages = c("data.table", "dplyr"),
    iteration = "group"
  ),
  tar_target(gwas_results_file,
    command = do_gwas(
      data = sample_sheet_qc, # phenotypes
      model = gwas_models,
      vcfs = ga_imputed_vcf_crossmap, # VCFs from tar_crossmap
      vep = veb_symbol_file, # VEP annotation from tar_vep
      path = file.path(output_directory, "gwas"),
      bin_path = list(
        bcftools = bcftools,
        plink2 = plink2
      )
    ),
    pattern = map(gwas_models),
    iteration = "list",
    packages = c("here", "data.table"),
    format = "file"
  )
)}

list(
  tar_setup,
  tar_sample_sheet_qc,
  tar_gwas
)