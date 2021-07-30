message(timestamp(quiet = TRUE))
### Project Setup ==================================================================================
library(here)
project_name <- sub("(.*)_*\\..*", "\\1", list.files(here(), pattern = ".Rproj$"))


### Input ==========================================================================================
epic <- here("outputs", "02-qc_idats")
epic_snps <- here("outputs", "02-qc_idats_snps")
omni <- here("outputs", "03-qc_plink")
ethnicity <- here("outputs", "05-ethnicity")
omni_hg38 <- here("outputs", "06-omni_to_hg38")


### Output =========================================================================================
output_directory <- normalizePath(file.path("/disks/DATA/Projects", project_name, "QC"))
sapply(
  X = file.path(output_directory, c("EPIC", "Omni2.5", "RNA-seq")),
  FUN = dir.create, recursive = TRUE, showWarnings = FALSE, mode = "0775"
)


### Save EPIC array ================================================================================
for (ifile in c("EPIC_QC_betavalues.csv.gz", "EPIC_QC_mset.rds", "EPIC_QC_phenotypes.csv")) {
  file.copy(
    from = file.path(epic, ifile),
    to = file.path(output_directory, "EPIC", ifile),
    overwrite = TRUE
  )
}
file.copy(
  from = here("reports", paste0(project_name, "_EPIC_QC.html")),
  to = file.path(output_directory, "EPIC", paste0(project_name, "_EPIC_QC.html")),
  overwrite = TRUE
)

for (ifile in c("EPIC_QC_betavalues.csv.gz", "EPIC_QC_mset.rds", "EPIC_QC_phenotypes.csv")) {
  file.copy(
    from = file.path(epic_snps, ifile),
    to = file.path(output_directory, "EPIC", gsub("^([^.]*)\\.(.*)$", "\\1_SNPs.\\2", ifile)),
    overwrite = TRUE
  )
}
file.copy(
  from = here("reports", paste0(project_name, "_EPIC_QC_SNPs.html")),
  to = file.path(output_directory, "EPIC", paste0(project_name, "_EPIC_QC_SNPs.html")),
  overwrite = TRUE
)


### Save Omni2.5 array =============================================================================
for (ifile in c(
  file.path("plink_qc", paste0(project_name, "_QC.bed")),
  file.path("plink_qc", paste0(project_name, "_QC.bim")),
  file.path("plink_qc", paste0(project_name, "_QC.fam")),
  file.path("vcf_qc", paste0(project_name, "_QC.vcf.gz")),
  file.path("vcf_qc", paste0(project_name, "_QC.vcf.gz.csi"))
)) {
  file.copy(
    from = file.path(omni, , ifile),
    to = file.path(output_directory, "Omni2.5", ifile),
    overwrite = TRUE
  )
}
file.copy(
  from = here("reports", paste0(project_name, "_OMNI2.5_QC.html")),
  to = file.path(output_directory, "Omni2.5", paste0(project_name, "_OMNI2.5_QC.html")),
  overwrite = TRUE
)
file.copy(
  from = file.path(omni, , "QC_exclusion.xlsx"),
  to = file.path(output_directory, "Omni2.5", paste0(project_name, "_QC_exclusion.xlsx")),
  overwrite = TRUE
)


### Imputed files ==================================================================================
if (file.exists(file.path(output_directory, "Omni2.5", "vcf_imputed"))) {
  unlink(
    list.files(
      path = file.path(output_directory, "Omni2.5", "vcf_imputed"), 
      pattern = "csi$|tbi$", 
      full.names = TRUE
    )
  )
  sapply(
    X = list.files(
      path = file.path(output_directory, "Omni2.5", "vcf_imputed"), 
      pattern = "vcf.gz$", 
      full.names = TRUE
    ),
    FUN = function(ivcf) system(paste("tabix -f -p vcf", ivcf))
  )
  identical(
    length(
      list.files(
        path = file.path(output_directory, "Omni2.5", "vcf_imputed"), 
        pattern = "vcf.gz$", 
        full.names = TRUE
      )
    ),
    length(
      list.files(
        path = file.path(output_directory, "Omni2.5", "vcf_imputed"), 
        pattern = "tbi$", 
        full.names = TRUE
      )
    )
  )
}
file.copy(
  from = here("reports", paste0(project_name, "_OMNI2.5_IMPUTE_QC.html")),
  to = file.path(output_directory, "Omni2.5", paste0(project_name, "_OMNI2.5_IMPUTE_QC.html")),
  overwrite = TRUE
)


### Ethnicity ======================================================================================
file.copy(
  from = file.path(ethnicity, , paste0(project_name, "_ethnicity.csv")),
  to = file.path(output_directory, "Omni2.5", paste0(project_name, "_ethnicity.csv")),
  overwrite = TRUE
)
file.copy(
  from = file.path(ethnicity, , paste0(project_name, "_ethnicity.pdf")),
  to = file.path(output_directory, "Omni2.5", paste0(project_name, "_ethnicity.pdf")),
  overwrite = TRUE
)


### Imputed files GRCh38 ===========================================================================
file.copy(
  from = list.files(file.path(omni_hg38), full.names = TRUE),
  to = file.path(output_directory, "Omni2.5", "vcf_imputed_hg38_ucsc"),
  overwrite = TRUE
)
if (file.exists(file.path(output_directory, "Omni2.5", "vcf_imputed_hg38_ucsc"))) {
  unlink(
    list.files(
      path = file.path(output_directory, "Omni2.5", "vcf_imputed_hg38_ucsc"), 
      pattern = "csi$|tbi$", 
      full.names = TRUE
    )
  )
  sapply(
    X = list.files(
      path = file.path(output_directory, "Omni2.5", "vcf_imputed_hg38_ucsc"), 
      pattern = "vcf.gz$", 
      full.names = TRUE
    ),
    FUN = function(ivcf) system(paste("tabix -f -p vcf", ivcf))
  )
  identical(
    length(
      list.files(
        path = file.path(output_directory, "Omni2.5", "vcf_imputed_hg38_ucsc"), 
        pattern = "vcf.gz$", 
        full.names = TRUE
      )
    ),
    length(
      list.files(
        path = file.path(output_directory, "Omni2.5", "vcf_imputed_hg38_ucsc"), 
        pattern = "tbi$", 
        full.names = TRUE
      )
    )
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
  list.files(output_directory, full.names = TRUE, recursive = TRUE, all.files = TRUE), 
  mode = "0555", use_umask = FALSE
)


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
message(timestamp(quiet = TRUE))
