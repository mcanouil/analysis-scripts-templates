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
  tar_target(genome_assembly, "GRCh38"),
  tar_target(ensembl_version, "104"),
  tar_target(ensembl_species, "homo_sapiens"),
  tar_target(vep_cache,
    command = c(
      "server" = "/media/Data/ExternalData/vep_data", 
      "docker" = "/disks/DATA/ExternalData/vep_data"
    )
  )
)}


### targets ========================================================================================
tar_vep <- {list(
  tar_target(snps_locations, "snps_locations.txt.gz", format = "file"), # CHR, START, END, A/C, + tab delimited
  tar_target(vep_symbol, 
    command = get_symbol_vep(
      input = snps_locations,
      output_directory = output_directory,
      genome_assembly = genome_assembly,
      ensembl_version = ensembl_version,
      ensembl_species = ensembl_species,
      vep_cache = vep_cache
    ),
    packages = c("here", "data.table")
  ),
  tar_target(veb_symbol_file,
    command = format_symbol_vep(vep_symbol),
    format = "file"
  )
)}

list(
  tar_setup,
  tar_vep
)
