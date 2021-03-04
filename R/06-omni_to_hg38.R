### Project Setup ==================================================================================
library(here)
project_name <- gsub("(.*)_.*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "06-omni_to_hg38")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

vcf_directory <- file.path("/disks/DATA/ExternalData", project_name, "QC/Omni2.5/vcf_imputed")


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(dgapaq)
})


### Analysis =======================================================================================
convert_assembly(
  input_directory = vcf_directory,
  output_directory = output_directory,
  ref_fasta = "/disks/DATA/ExternalData/1kg/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
  chain_file = "/disks/DATA/ExternalData/1kg/chain_files/hg19ToHg38.over.chain.gz",
  bin_path = list(
    crossmap = "/usr/local/bin/CrossMap.py",
    bcftools = "/usr/bin/bcftools",
    tabix = "/usr/bin/tabix",
    bgzip = "/usr/bin/bgzip"
  ),
  nb_cores = 22
)


### Archive ========================================================================================
# if (!Sys.getenv("USER") %in% c("root", "") && file.exists("~/.fex/id")) {
#   local({
#     owd <- getwd()
#     setwd(normalizePath(output_directory))
#     archive_name <- file.path(
#       normalizePath(output_directory),
#       paste0(
#         format(Sys.Date(), format = "%Y%m%d"), "_",
#         project_name, "_", 
#         gsub("[0-9]+\\-", "", basename(output_directory)), ".zip"
#       )
#     )
#     zip(archive_name, files = list.files())
#     fex_out <- system(paste("fexsend", archive_name, "."), intern = TRUE)
#     unlink(archive_name)
#     setwd(owd)
#   })
# }


### Set chmod ======================================================================================
Sys.chmod(
  list.files(output_directory, full.names = TRUE),
  mode = "0775", use_umask = FALSE
)
Sys.chmod(
  list.files(output_directory, full.names = TRUE, recursive = TRUE, all.files = TRUE),
  mode = "0775", use_umask = FALSE
)
invisible(system(paste("chgrp -R staff", output_directory), intern = TRUE))


### Complete =======================================================================================
message("Success!", appendLF = TRUE)

