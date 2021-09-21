message(timestamp(quiet = TRUE))
### Project Setup ==================================================================================
library(here)
project_name <- sub("(.*)_*\\..*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "14-mqtl")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

invisible(sapply(c("epic", "covariates", "fastqtl", "fastqtl_annotated", "fastqtl_combined"), function(x) {
  assign(x = paste0("output_", x), value = file.path(output_directory, x), envir = .GlobalEnv)
  dir.create(file.path(output_directory, x), recursive = TRUE, showWarnings = FALSE, mode = "0775")
}))

vep_directory <- here("outputs", "08-vep_vcf_docker")

data_directory <- file.path("/disks/DATA/Projects", project_name, "QC")

phenotype <- here("docs", "phenotype.xlsx")
ethnicity <- file.path(data_directory, "Omni2.5", paste0(project_name, "_ethnicity.csv"))
exclusion <- file.path(data_directory, "Omni2.5", paste0(project_name, "_QC_exclusion.xlsx"))


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(data.table)
  library(future)
  library(future.apply)
  library(future.callr)
})


### project setup ==================================================================================
plan(future.callr::callr, workers = 11)
message(sprintf("Number of workers: %d", future::nbrOfWorkers()))


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

epic_phenotype <- fread(
  file = file.path(data_directory, "EPIC", "EPIC_QC_phenotypes.csv"),
  colClasses = c("Sample_ID" = "character")
)
sample_sheet_qc <- merge(
  x = sample_sheet_qc,
  y = epic_phenotype[j = .SD, .SDcols = grep("^CellT|^Sample_ID|^Sentrix_", names(epic_phenotype))],
  by.x = "#IID",
  by.y = "Sample_ID"
)

relatedness <- setDT(read_excel(exclusion, 2))[gsub("_.*", "", IID1) != gsub("_.*", "", IID2)]
relatedness[j = paste0("vcf_id", 1:2) := list(paste0(FID1, "_", IID1), paste0(FID2, "_", IID2))]

sample_sheet_qc <- sample_sheet_qc[ # not related
  !vcf_id %in% unique(unlist(relatedness[j = c("vcf_id1", "vcf_id2")], use.names = FALSE))
][ # European descent
  super_pop_closest %in% "EUR"
]


## Omni2.5 data ====================================================================================
vcfs <- list.files(
  path = file.path(vep_directory, "vcfs_qc"),
  pattern = "[^X].vcf.gz$",
  full.names = TRUE
)
names(vcfs) <- gsub("_qc.vcf.gz$", "", basename(vcfs))


### EPIC data ======================================================================================
beta_matrix <- fread(
  file = file.path(data_directory, "EPIC", "EPIC_QC_betavalues_SNPs.csv.gz"),
  header = TRUE
)
beta_matrix <- (function(x) log2(x) - log2(1 - x))(as.matrix(beta_matrix, "cpg_id"))[, sample_sheet_qc[["Sample_ID"]]]

epic_qc_annot <- Locations
epic_qc_annot <- epic_qc_annot[intersect(rownames(beta_matrix), rownames(epic_qc_annot)), ]
epic_qc_annot <- na.exclude(as.data.frame(epic_qc_annot))
epic_qc_annot <- epic_qc_annot[epic_qc_annot[["chr"]] %in% sprintf("chr%d", 1:22), ]
epic_qc_annot[["chr"]] <- as.numeric(gsub("chr", "", epic_qc_annot[["chr"]]))
epic_qc <- beta_matrix[rownames(epic_qc_annot), as.character(sample_sheet_qc[["Sample_ID"]])]
colnames(epic_qc) <- sample_sheet_qc[["vcf_id"]]
epic_qc <- merge(
  x = as.data.table(epic_qc, keep.rownames = "CpG"),
  y = as.data.table(epic_qc_annot, keep.rownames = "CpG"),
  by = "CpG"
)
epic_qc[j = (c("#Chr", "start", "end", "grp")) := list(chr, pos, pos, chr)]
epic_qc <- epic_qc[j = .SD, .SDcols = c("grp", "#Chr", "start", "end", "CpG", sample_sheet_qc[["vcf_id"]])]
epic_qc <- epic_qc[order(`#Chr`, start)]

epic_qc[
  j = (function(x, y) {
    fwrite(x, file = file.path(output_epic, sprintf("chr%02d.bed", unique(y))), quote = FALSE, sep = "\t")
    system(paste("bgzip -f", file.path(output_epic, sprintf("chr%02d.bed", unique(y)))))
    system(paste("tabix -p bed -f", file.path(output_epic, sprintf("chr%02d.bed.gz", unique(y)))))
    return(TRUE)
  })(.SD, `#Chr`),
  by = grp
]


### mQTL analysis (fastQTL) =======================================================================
file_con <- gzfile(file.path(tempdir(), "nominal_header.txt.gz"), "w")
cat(
  c("cpg_id", "rs_id", "distance_bp", "pvalue", "slope\n"),
  sep = " ",
  file = file_con
)
close(file_con)

file_con <- gzfile(file.path(tempdir(), "permutation_header.txt.gz"))
cat(
  c(
    "cpg_id", "variants_cis", "mle_shape1_beta", "mle_shape2_beta",
    "dummy", "best_rs_id", "distance_bp", "pvalue", "slope",
    "permutation_pvalue", "downstream_pvalue\n"
  ),
  sep = " ",
  file = file_con
)
close(file_con)

for (ichr in sprintf("chr%02d", 1:22)) {
  local({
    for (ianalysis in c("nominal", "permutation")) {
      n_chunk <- 20
      cis_window <- 500000 # 500 Kb
      future_lapply(
        X = 1:n_chunk,
        future.globals = FALSE,
        vcfs = vcfs,
        output_epic = output_epic,
        output_covariates = output_covariates,
        output_fastqtl = output_fastqtl,
        cis_window = cis_window,
        ianalysis = ianalysis,
        ichr = ichr,
        n_chunk = n_chunk
        cis_window = cis_window,
        FUN = function(ichunk, vcfs, output_epic, output_covariates, output_fastqtl, cis_window, ianalysis, ichr, n_chunk) {
          system(paste("fastQTL",
            "--silent",
            "--seed", 20200331,
            "--vcf", vcfs[ichr],
            "--bed", file.path(output_epic, paste0(ichr, ".bed.gz")),
            "--cov", file.path(output_covariates, "covariates.txt.gz"),
            "--window", cis_window,
            ifelse(ianalysis == "nominal", "", "--permute 1000 10000"),
            "--chunk", ichunk, n_chunk,
            "--out", file.path(output_fastqtl, sprintf("%s_%s_%03d.txt.gz", ianalysis, ichr, ichunk))
          ))
        }
      )

      system(paste(
        "zcat",
        file.path(tempdir(), sprintf("%s_header.txt.gz", ianalysis)),
        file.path(output_fastqtl, sprintf("%s_%s_*.txt.gz", ianalysis, ichr)),
        "| bgzip -c >", file.path(output_fastqtl, sprintf("%s_%s_%s.txt.gz", project_name, ianalysis, ichr))
      ))

      unlink(list.files(
        path = output_fastqtl,
        pattern = sprintf("%s_%s_[0-9]*.txt.gz", ianalysis, ichr),
        full.names = TRUE
      ))

      fwrite(
        x = Reduce(
          f = function(x, y) merge(x, y, by = "cpg_id", all.x = TRUE),
          x = list(
            fread(file.path(output_fastqtl, sprintf("%s_%s_%s.txt.gz", project_name, ianalysis, ichr))),
            as.data.table(epic_qc_annot, keep.rownames = "cpg_id"),
            as.data.table(Other, keep.rownames = "cpg_id")[
              j = list(
                UCSC_RefGene_Name = paste(unique(tstrsplit(UCSC_RefGene_Name, split = ";")), collapse = ";")
              ),
              by = "cpg_id"
            ]
          )
        ),
        file = file.path(output_fastqtl_annotated, sprintf("%s_%s_%s.txt.gz", project_name, ianalysis, ichr))
      )
    }
  })
}


### Combine Results ================================================================================
for (ianalysis in c("nominal", "permutation")) {
  local({
    fwrite(
      x = rbindlist(future_lapply(
        X = sprintf("chr%02d", 1:22),
        future.globals = FALSE,
        future.packages = "data.table",
        output_fastqtl_annotated = output_fastqtl_annotated,
        project_name = project_name,
        ianalysis = ianalysis,
        FUN = function(ichr, output_fastqtl_annotated, project_name, ianalysis) {
          fread(file.path(output_fastqtl_annotated, sprintf("%s_%s_%s.txt.gz", project_name, ianalysis, ichr)))
        }
      )),
      file = file.path(output_fastqtl_combined, sprintf("%s_%s.txt.gz", project_name, ianalysis))
    )
  })
}


### Clean Directory ================================================================================
unlink(c(
  gzfile(file.path(tempdir(), "nominal_header.txt.gz")),
  gzfile(file.path(tempdir(), "permutation_header.txt.gz"))
))
unlink(setdiff(
  list.files(output_directory, full.names = TRUE),
  list.files(output_directory, pattern = "fastqtl_combined", full.names = TRUE)
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
message(timestamp(quiet = TRUE))
