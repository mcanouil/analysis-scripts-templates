### global libraries ===============================================================================
library(targets)
library(tarchetypes)
library(here)
library(data.table)
# library(future)
# library(future.callr)

# targets::tar_renv(extras = "visNetwork", path = "scripts/_dependencies.R")


### project setup ==================================================================================
# Functions/scripts required: tar-crossmap.R
invisible(sapply(list.files(here("scripts"), pattern = "^tar-.*R$", full.names = TRUE), source, echo = FALSE))

# plan(future.callr::callr, workers = 40)
# plan(multicore, workers = 40)
# message(sprintf("Number of workers: %d", future::nbrOfWorkers()))
# setDTthreads(threads = 1)


### targets setup ==================================================================================
tar_setup <- list(
  tar_target(crossmap, "/usr/local/bin/CrossMap.py", format = "file"),
  tar_target(bcftools, "/usr/bin/bcftools", format = "file"),
  tar_target(bin_path, list(bcftools = bcftools, crossmap = crossmap)),
  tar_target(reference_fasta,
    command = "/disks/DATA/ExternalData/1kg/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
    format = "file"
  ),
  tar_target(chain_file,
    command = "/disks/DATA/ExternalData/1kg/chain_files/hg19ToHg38.over.chain.gz",
    format = "file"
  )
)


### targets ========================================================================================
tar_crossmap <- {list(
  tar_target(ga_imputed_vcf_crossmap,
    command = {
      if (ga_save_qc) {
        do_crossmap(
          path = ga_imputed_vcf, # path to directory or chr VCF files
          output_directory = file.path(ga_export_directory, "vcf_imputed_grch38"),
          ref_fasta = reference_fasta,
          chain_file = chain_file,
          bin_path = bin_path
        )
      }
    },
    packages = "future.apply",
    format = "file"
  )
)}

list(
  tar_setup,
  tar_crossmap
)
