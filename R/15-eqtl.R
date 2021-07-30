message(timestamp(quiet = TRUE))
### Project Setup ==================================================================================
library(here)
project_name <- sub("(.*)_*\\..*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "15-eqtl")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

vep_directory <- here("outputs", "08-vep_vcf_docker")

run_directory <- "/disks/RUN/Run_XXX/Output/RSEM"

phenotype <- here("docs", "phenotype.xlsx")
ethnicity <- file.path(data_directory, "Omni2.5", paste0(project_name, "_ethnicity.csv"))
exclusion <- file.path(data_directory, "Omni2.5", paste0(project_name, "_QC_exclusion.xlsx"))

seed <- as.numeric(Sys.Date())
cis_window <- 500


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(parallel)
  library(tximport)
  library(DESeq2)
  library(matrixStats)
  library(glue)
  library(biomaRt)
  library(httr)
  library(data.table)
})


### Setup biomaRt ==================================================================================
set_config(config(ssl_verifypeer = FALSE)) # Fix SSL check error

species <- "hsapiens_gene_ensembl"
build <- 38
version <- 103

get_mart <- quote(useEnsembl(
  biomart = "ensembl", 
  dataset = species, 
  version = version, 
  GRCh = if (build == 37) build else NULL
))

mart <- try(eval(get_mart), silent = TRUE)
if (inherits(mart, "try-error")) mart <- eval(get_mart)
ensembl_version <- sprintf("GRCh%d-%d", build, version)


### Analysis =======================================================================================
for (rna_level in c("genes", "isoforms")) { 
  local({
    rna_level_name <- unname(c("genes" = "gene", "isoforms" = "transcript")[rna_level])
    
    output_directory <- here("outputs", "14-eqtl", rna_level)
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")
    
    invisible(sapply(c("rnaseq", "covariates", "fastqtl", "fastqtl_annotated", "fastqtl_combined"), function(x) {
      assign(x = paste0("output_", x), value = file.path(output_directory, x), envir = .GlobalEnv)
      dir.create(file.path(output_directory, x), recursive = TRUE, showWarnings = FALSE, mode = "0775")
    }))
    
    
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
      "bmi" = bmi,
      "rnaseq_path" = file.path(run_directory, "...")
    )]
    
    relatedness <- setDT(read_excel(exclusion, 2))[gsub("_.*", "", IID1) != gsub("_.*", "", IID2)]
    relatedness[j = paste0("vcf_id", 1:2) := list(paste0(FID1, "_", IID1), paste0(FID2, "_", IID2))]
    
    sample_sheet_qc <- sample_sheet_qc[ # not related
      !vcf_id %in% unique(unlist(relatedness[j = c("vcf_id1", "vcf_id2")], use.names = FALSE))
    ][ # European descent
      super_pop_closest %in% "EUR"
    ]
    
    fwrite(
      x = transpose(sample_sheet_qc[j = c("vcf_id", "sex", "age", "bmi", "t2d", "PC01", "PC02")], keep.names = "row"),
      file = file.path(output_covariates, "covariates.txt.gz"),
      quote = FALSE, col.names = FALSE, sep = "\t"
    )
    
    
### Omni2.5 data ===================================================================================
    vcfs <- list.files(
      path = file.path(vep_directory, "vcfs_qc"),
      pattern = "[^X].vcf.gz$",
      full.names = TRUE
    )
    names(vcfs) <- gsub("_qc.vcf.gz$", "", basename(vcfs))
    
    
### RNAseq data ====================================================================================
    rsem_files <- setNames(sample_sheet_qc[["rsem_file"]], sample_sheet_qc[["vcf_id"]])
    txi_counts <- tximport(
      files = rsem_files, 
      type = "rsem", 
      txIn = rna_level == "isoforms", 
      txOut = rna_level == "isoforms", 
      countsFromAbundance = "no"
    )
    txi_counts$length[txi_counts$length == 0] <- 1
    txi_counts[["counts"]] <- round(txi_counts[["counts"]], digits = 0)
    class(txi_counts[["counts"]]) <- "integer"
    gene_matrix <- txi_counts[["counts"]][
      rowVars(txi_counts[["counts"]]) != 0 & 
        rowMeans(txi_counts[["counts"]]) > 1 & 
        rowMedians(txi_counts[["counts"]]) > 0, 
    ]
    gene_matrix <- as.data.table(vst(gene_matrix))[j = ensembl_id := rownames(gene_matrix)]
    
    if (file.exists(file.path(output_rnaseq, "rnaseq_qc_map.csv.gz"))) {
      rnaseq_qc_map <- fread(file.path(output_rnaseq, "rnaseq_qc_map.csv.gz"))
    } else {
      rnaseq_qc_map <- setDT(getBM(
        attributes = c(glue("ensembl_{rna_level_name}_id"), "chromosome_name", "start_position", "end_position", "external_gene_name"),
        filters = glue("ensembl_{rna_level_name}_id"),
        values = list(gene_matrix[["ensembl_id"]]),
        mart = mart
      ))[
        j = lapply(.SD, function(x) {
          out <- paste(setdiff(unique(x), ""), collapse = ";")
          fifelse(out == "", NA_character_, out)
        }),
        by = c(glue("ensembl_{rna_level_name}_id"))
      ][j = ensembl_version := ensembl_version]
      fwrite(rnaseq_qc_map, file = file.path(output_rnaseq, "rnaseq_qc_map.csv.gz"))
    }
    
    rnaseq_qc <- merge(
      x = gene_matrix, 
      y = rnaseq_qc_map[j = -"external_gene_name"], 
      by.x = "ensembl_id",
      by.y = glue("ensembl_{rna_level_name}_id")
    )[
      chromosome_name %in% 1:22
    ][
      j = c("grp", "#Chr", "start", "end") := list(
        as.numeric(chromosome_name),
        as.numeric(chromosome_name),
        start_position,
        end_position
      )
    ][
      j = .SD,
      .SDcols = c("grp", "#Chr", "start", "end", "ensembl_id", setdiff(colnames(gene_matrix), "ensembl_id"))
    ][order(grp, start)]
    
    rnaseq_qc[
      j = (function(x, y) {
        fwrite(x, file = file.path(output_rnaseq, sprintf("chr%02d.bed", unique(y))), quote = FALSE, sep = "\t")
        system(paste("bgzip -f", file.path(output_rnaseq, sprintf("chr%02d.bed", unique(y)))))
        system(paste("tabix -p bed -f", file.path(output_rnaseq, sprintf("chr%02d.bed.gz", unique(y)))))
        return(TRUE)
      })(.SD, `#Chr`), 
      by = grp
    ]
    
    
### eQTL analysis (fastQTL) ========================================================================
    file_con <- gzfile(file.path(tempdir(), "nominal_header.txt.gz"), "w")
    cat(
      c("ensembl_id", "rs_id", "distance_bp", "pvalue", "slope\n"),
      sep = " ",
      file = file_con
    )
    close(file_con)
    
    file_con <- gzfile(file.path(tempdir(), "permutation_header.txt.gz"))
    cat(
      c(
        "ensembl_id", "variants_cis", "mle_shape1_beta", "mle_shape2_beta",
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
          mclapply(1:n_chunk, mc.cores = n_chunk, function(ichunk) {
            system(paste("fastQTL",
              "--silent",
              "--seed", seed,
              "--vcf", vcfs[ichr],
              "--bed", file.path(output_rnaseq, paste0(ichr, ".bed.gz")),
              "--cov", file.path(output_covariates, "covariates.txt.gz"),
              "--window", cis_window,
              ifelse(ianalysis == "nominal", "", "--permute 1000 10000"),
              "--chunk", ichunk, n_chunk,
              "--out", file.path(output_fastqtl, sprintf("%s_%s_%03d.txt.gz", ianalysis, ichr, ichunk))
            ))
          })
      
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
            x = merge(
              x = fread(file.path(output_fastqtl, sprintf("%s_%s_%s.txt.gz", project_name, ianalysis, ichr))),
              y = rnaseq_qc_map,
              by.x = "ensembl_id",
              by.y = glue("ensembl_{rna_level_name}_id"),
              all.x = TRUE
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
          x = rbindlist(mclapply(
            X = sprintf("chr%02d", 1:22),
            mc.cores = 11, 
            mc.preschedule = FALSE,
            FUN = function(ichr) {
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
  })
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
