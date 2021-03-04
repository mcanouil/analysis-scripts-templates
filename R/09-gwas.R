### Project Setup ==================================================================================
library(here)
project_name <- gsub("(.*)_.*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "09-gwas")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

invisible(sapply(c("lm", "lmqc", "annot", "lmqc_annot"), function(x) {
  assign(x = paste0("output_", x), value = file.path(output_directory, x), envir = .GlobalEnv)
  dir.create(file.path(output_directory, x), recursive = TRUE, showWarnings = FALSE, mode = "0775")
}))

phenotype <- here("docs", "phenotype.xlsx")
vep_directory <- here("outputs", "08-vep_vcf_docker")

data_directory <- file.path("/disks/DATA/Projects", project_name, "QC")
ethnicity <- file.path(data_directory, "Omni2.5", paste0(project_name, "_ethnicity.csv"))
exclusion <- file.path(data_directory, "Omni2.5", paste0(project_name, "_QC_exclusion.xlsx"))

plink_url <- "http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip"


### Load packages ==================================================================================
suppressPackageStartupMessages({
  library(glue)
  library(parallel)
  library(data.table)
  library(readxl)
})


### Get PLINK2 =====================================================================================
if (!file.exists(file.path(output_directory, "plink2"))) {
  download.file(
    url = plink_url,
    destfile = file.path(output_directory, "plink2.zip")
  )
  unzip(file.path(output_directory, "plink2.zip"), exdir = output_directory)
  unlink(file.path(output_directory, "plink2.zip"))
  Sys.chmod(file.path(output_directory, "plink2"), "0777")
}


### Get QCed Phenotypes ============================================================================
sample_sheet_qc <- merge(
  x = setDT(read_excel(phenotype)),
  y = fread(ethnicity),
  by = "iid"
)[j = `:=`(
  "Sample_ID" = iid,
  "#IID" = iid,
  "vcf_id" = iid,
  "group" = as.integer(groupe) + 1L,
  "sex" = c("M" = 1L, "F" = 2L)[sex],
  "age" = age,
  "bmi" = bmi
)]

relatedness <- setDT(read_excel(exclusion, 2))[gsub("_.*", "", IID1) != gsub("_.*", "", IID2)]
relatedness[j = paste0("vcf_id", 1:2) := list(paste0(FID1, "_", IID1), paste0(FID2, "_", IID2))]

sample_sheet_qc <- sample_sheet_qc[ # not related
  !vcf_id %in% unique(unlist(relatedness[j = c("vcf_id1", "vcf_id2")], use.names = FALSE))
][ # European descent
  super_pop_closest %in% "EUR"
]


### Traits to be analysed ==========================================================================
traits <- c(
  "Control/Case" = "group"
)

covariates <- list(
  "group" = c("sex", "age", "bmi", sprintf("PC%02d", 1:2))
)


### Format VCF =====================================================================================
#' Criteria
#' + Samples
#'   * European (1k Genome)
#'   * Call rate greater than 0.90
#'   * Sex consistency
#'   * Heterozygosity rate lower than four times the SD from the mean
#' + Variants
#'   * INFO greater than 0.8 (a posteriori)
#'   * SNPs
#'   * Hardy-Weinberg equilibrium test above 0.005 (a posteriori)
#'   * maf greater than 0.05 (a posteriori)
cat(
  sample_sheet_qc[["vcf_id"]],
  file = file.path(output_directory, "samples_to_keep.txt"),
  sep = "\n"
)
vcfs <- list.files(
  path = file.path(vep_directory, "vcfs_qc"),
  pattern = "[^X].vcf.gz$",
  full.names = TRUE
)
names(vcfs) <- gsub("_qc.vcf.gz$", "", basename(vcfs))

if (length(list.files(output_annot)) != 4 * 22) {
  local({
    mclapply(
      X = names(vcfs),
      mc.cores = 11,
      mc.preschedule = FALSE,
      FUN = function(ivcf) {
        sapply(
          X = paste(
            "vcftools --gzvcf", vcfs[ivcf],
            c("--get-INFO 'INFO'", "--freq", "--hardy", "--missing-site"),
            "--out", file.path(output_annot, ivcf)
          ),
          FUN = system
        )
        system(paste(
          'sed -i s/"{ALLELE:FREQ}"/frq1"\t"frq2/',
          file.path(output_annot, paste0(ivcf, ".frq"))
        ))
      }
    )
  })
}


### Analyses =======================================================================================
for (trait in traits) {
  local({
    dir.create(file.path(output_lm, trait), recursive = TRUE, showWarnings = FALSE, mode = "0775")
    
    keep_samples <- sample_sheet_qc[!is.na(get(trait)), c("#IID")]
    fwrite(
      x = keep_samples, 
      file = file.path(output_lm, trait, glue("{trait}.samples")),
      sep = "\t"
    )
    
    fwrite(
      x = sample_sheet_qc[`#IID` %in% keep_samples[[1]], .SD, .SDcols = c("#IID", trait)], 
      file = file.path(output_lm, trait, glue("{trait}.pheno")),
      sep = " "
    )

    fwrite(
      x = sample_sheet_qc[`#IID` %in% keep_samples[[1]], .SD, .SDcols = c("#IID", setdiff(covariates[[trait]], "sex"))],
      file = file.path(output_lm, trait, glue("{trait}.cov")),
      sep = " "
    )
  
    fwrite(
      x = sample_sheet_qc[`#IID` %in% keep_samples[[1]], .(`#IID`, SEX = sex)],
      file = file.path(output_lm, trait, glue("{trait}.sex")),
      sep = " "
    )
    
    model <- if (length(unique(sample_sheet_qc[`#IID` %in% keep_samples[[1]]][[trait]])) == 2) {
      "logistic"
    } else {
      "linear"
    }
    
    for (ivcf in names(vcfs)) {
      local({
        output_file <- file.path(
          output_lm,
          trait,
          paste(model, ivcf, sep = "_")
        )
        system(paste(
          file.path(output_directory, "plink2"),
          "--vcf", vcfs[ivcf], "dosage=DS",
          "--mach-r2-filter",
          "--threads 120",
          "--glm sex",
          "--keep", file.path(output_lm, trait, glue("{trait}.samples")),
          "--update-sex", file.path(output_lm, trait, glue("{trait}.sex")),
          "--pheno", file.path(output_lm, trait, glue("{trait}.pheno")),
          "--covar", file.path(output_lm, trait, glue("{trait}.cov")),
          "--covar-variance-standardize",
          "--out", output_file
        ))
        system(paste("bgzip --force", paste0(output_file, glue(".{trait}.glm.{model}"))))
      })
    }
  })
}


### QC  ============================================================================================
for (trait in traits) {
  local({
    ntrait <- names(traits)[traits %in% trait]
    ifile <- any(grepl("linear_", list.files(file.path(output_lm, trait))))
    model <- if (ifile) "linear" else "logistic"

    res_gwas <- rbindlist(
      lapply(
        X = names(vcfs),
        FUN = function(ivcf) {
          plink_res <- setnames(
            x = fread(
              file = file.path(output_lm, trait, paste0(model, "_", ivcf, ".", trait, ".glm.", model, ".gz")), 
              header = TRUE
            ),
            old = "#CHROM",
            new = "CHR",
            skip_absent = TRUE
          )
    
          plink_annot <- Reduce(
            f = function(x, y) merge(x, y, by = c("CHR", "POS")),
            x = lapply(
              X = list.files(path = output_annot, pattern = ivcf, full.names = TRUE),
              FUN = function(.x) setnames(fread(.x), "CHROM", "CHR", skip = TRUE)
            )
          )
    
          trait_res <- merge(
            x = plink_res,
            y = plink_annot,
            by = c("CHR", "POS", "REF", "ALT"),
            all.x = TRUE
          )[, 
            (c("REF_FRQ", "ALT_FRQ")) := 
              list(as.numeric(gsub(".:", "", frq1)), as.numeric(gsub(".:", "", frq2)))
          ][
            is.finite(ChiSq_HWE) & 
              !is.na(P) & 
              TEST == "ADD" & 
              INFO >= 0.8 & 
              REF_FRQ >= 0.05 & 
              ALT_FRQ >= 0.05 & 
              P_HWE >= 0.005,
            -c("N_DATA", "frq1", "frq2", "OBS(HOM1/HET/HOM2)", "E(HOM1/HET/HOM2)", "P_HET_DEFICIT", "P_HET_EXCESS")
          ]
          
          multiallele_snp <- trait_res[duplicated(trait_res, by = c("CHR", "POS", "ID")), ID]
   
          unique(trait_res[!ID %in% multiallele_snp])[, FDR := p.adjust(P, method = "BH")][]
        }
      )
    )
    fwrite(
      x = res_gwas,
      file = file.path(output_lmqc, glue("{project_name}_GWAS_{trait}_{model}.csv.gz"))
    )
  })
}


### Annotate with Symbol ===========================================================================
vep_annotation <- fread(file.path(vep_directory, "snps_vep_v102.0_GRCh38_formated.txt.gz"))

for (trait in traits) {
  local({
    ntrait <- names(traits)[traits %in% trait]
    ifile <- any(grepl("linear", list.files(path = output_lmqc, pattern = paste0("_", trait, "_"))))
    model <- if (ifile) "linear" else "logistic"
    
    res_gwas <- fread(file.path(output_lmqc, glue("{project_name}_GWAS_{trait}_{model}.csv.gz")))
    res_gwas[, chr_pos_ref_alt := paste0(CHR, "_", POS, "_", REF, "/", ALT)]
    res_annot <- merge(
      x = res_gwas, 
      y = vep_annotation, 
      by = c("chr_pos_ref_alt", "CHR", "POS"), 
      all.x = TRUE
    )
    res_annot[(ID != rsid), rsid := fifelse(rsid == "-" | grepl("^rs", ID), ID, rsid)]
    res_annot[, c("ID", "rsid", "Gene", "Symbol") := list(
      rsid,
      NULL,
      fifelse(Gene == "-", NA_character_, Gene),
      fifelse(Symbol == "", NA_character_, Symbol)
    )]
    
    fwrite(
      x = res_annot,
      file = file.path(output_lmqc_annot, glue("{project_name}_GWAS_{trait}_{model}.csv.gz"))
    )
  })
}


### Clean directory ================================================================================
unlink(setdiff(
  list.files(output_directory, full.names = TRUE),
  list.files(output_directory, pattern = "lmqc", full.names = TRUE)
), recursive = TRUE, force = TRUE)


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

