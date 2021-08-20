message(timestamp(quiet = TRUE))
### Project Setup ==================================================================================
library(here)
project_name <- sub("(.*)_*\\..*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "05-ethnicity")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

vcf_directory <- file.path("/disks/DATA/ExternalData", project_name, "QC/Omni2.5/vcf_imputed")


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(rain)
})


### Analyses =======================================================================================
if (!file.exists(file.path(output_directory, "all.bim"))) {
  estimate_ethnicity(
    cohort_name = project_name,
    input_vcfs = vcf_directory,
    input_type = "array",
    output_directory = output_directory,
    ref1kg_vcfs = normalizePath("/disks/DATA/ExternalData/1kg/hg19/Genotypes/RawVCF"),
    ref1kg_population = normalizePath("/disks/DATA/ExternalData/1kg/samples_description/integrated_call_samples_v3.20130502.ALL.panel"),
    ref1kg_maf = 0.01,
    splitted_by_chr = TRUE,
    quality_tag = "INFO",
    quality_threshold = 0.9,
    n_cores = 11,
    bin_path = list(
      vcftools = normalizePath("/usr/bin/vcftools"),
      bcftools = normalizePath("/usr/bin/bcftools"),
      bgzip = normalizePath("/usr/bin/bgzip"),
      tabix = normalizePath("/usr/bin/tabix"),
      plink = normalizePath("/usr/bin/plink1.9")
    )
  )
} else {
  compute_pca(
  	cohort_name = project_name, 
  	input_plink = file.path(output_directory, "all"), 
  	output_directory = output_directory, 
  	ref1kg_population = normalizePath("/disks/DATA/ExternalData/1kg/samples_description/integrated_call_samples_v3.20130502.ALL.panel")
  )
}


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
message(timestamp(quiet = TRUE))
