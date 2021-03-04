### Project Setup ==================================================================================
library(here)
project_name <- gsub("(.*)_.*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "10-ewas")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

phenotype <- here("docs", "phenotype.xlsx")

data_directory <- file.path("/disks/DATA/Projects", project_name, "QC")
ethnicity <- file.path(data_directory, "Omni2.5", paste0(project_name, "_ethnicity.csv"))
exclusion <- file.path(data_directory, "Omni2.5", paste0(project_name, "_QC_exclusion.xlsx"))


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(limma)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(DMRcate)
  library(DMRcatedata)
  library(glue)
  library(data.table)
  library(readxl)
})


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


### Traits to be analysed ==========================================================================
traits <- c(
  "Control/Case" = "group"
)

covariates <- list(
  "group" = c("sex", "age", "bmi", sprintf("PC%02d", 1:2))
)

cells <- list(
  "group" = grep("^CellT_", names(sample_sheet_qc), value = TRUE)
)


### EPIC data ======================================================================================
beta_matrix <- fread(
  file = file.path(data_directory, "EPIC", "EPIC_QC_betavalues.csv.gz"), 
  header = TRUE
)[j = .SD, .SDcols = c("cpg_id", sample_sheet_qc[["epic"]])]
beta_matrix <- setDF(beta_matrix[, -"cpg_id"], rownames = beta_matrix[["cpg_id"]])
beta_matrix <- as.matrix((function(x) log2(x) - log2(1 - x))(beta_matrix))


### Analyses =======================================================================================
for (trait in traits) {
  dir.create(file.path(output_directory, trait), recursive = TRUE, showWarnings = FALSE, mode = "0755")
  local({
    form <- as.formula(paste0("~ ", trait, " + ", paste(covariates[[trait]], collapse = " + ")))

    for (imodel in c("", "_cell")) {
      if (imodel == "_cell") {
        form <- update.formula(form, as.formula(glue("~ . + { paste(cells[[trait]], collapse = ' + ') }")))
      }

      pheno <- na.exclude(sample_sheet_qc[, .SD, .SDcols = c("ABOS_ID", "epic", all.vars(form))])

      limma_fit1 <- lmFit(
        object = beta_matrix[, as.character(pheno[["epic"]])], 
        design = model.matrix(object = form, data = pheno)
      )

      if (is.factor(sample_sheet_qc[[trait]]) & nlevels(sample_sheet_qc[[trait]]) > 2) {
        limma_fit2 <- eBayes(limma_fit1)
        limma_top1 <- lapply(
          X = paste0(trait, levels(sample_sheet_qc[[trait]])[-1]),
          .fit = limma_fit2,
          .trait = trait,
          .ref = levels(sample_sheet_qc[[trait]])[1],
          FUN = function(.coef, .fit, .trait, .ref) {
            out <- as.data.table(topTable(
              fit = .fit,
              coef = .coef,
              number = nrow(beta_matrix),
              adjust.method = "BH",
              p.value = 1,
              sort.by = "none"
            ), keep.rownames = "CpG")[
              j = `:=`(
                "se" = sqrt(.fit[["s2.post"]]) * .fit[["stdev.unscaled"]][, .coef],
                "adj.P.Val" = NULL,
                "B" = NULL,
                "fdr" = p.adjust(P.Value, method = "BH"),
                "contrast" = paste0(.trait, ": ", gsub(.trait, "", .coef), " Vs. ", .ref, " (ref)")
              )
            ]
            out
          }
        )
        limma_top2 <- apply(
          X = combn(levels(sample_sheet_qc[[trait]])[-1], 2),
          MARGIN = 2,
          .fit = limma_fit1,
          .trait = trait,
          FUN = function(icol, .fit, .trait) {
            out <- paste0(.trait, icol)
            coef <- paste0(out[2], "-", out[1])
            contrasts_fit <- makeContrasts(contrasts = coef, levels = .fit$design)
            attr(contrasts_fit, "dimnames") <- lapply(
              X = attr(contrasts_fit, "dimnames"),
              FUN = gsub, pattern = "^Intercept$", replacement = "(Intercept)"
            )
            limma_fit1b <- contrasts.fit(fit = .fit, contrasts = contrasts_fit)
            limma_fit2b <- eBayes(limma_fit1b)
            limma_top <- as.data.table(topTable(
              fit = limma_fit2b,
              coef = coef,
              number = nrow(beta_matrix),
              adjust.method = "BH",
              p.value = 1,
              sort.by = "none"
            ), keep.rownames = "CpG")[
              j = `:=`(
                "se" = sqrt(limma_fit2b[["s2.post"]]) * limma_fit2b[["stdev.unscaled"]][, coef],
                "adj.P.Val" = NULL,
                "B" = NULL,
                "fdr" = p.adjust(P.Value, method = "BH"),
                "contrast" = paste0(.trait, ": ", icol[2], " Vs. ", icol[1], " (ref)")
              )
            ]
            limma_top
          }
        )
        limma_top <- rbindlist(c(limma_top1, limma_top2))
      } else {
        limma_fit2 <- eBayes(limma_fit1)
        if (is.factor(sample_sheet_qc[[trait]])) {
          .levels <- levels(sample_sheet_qc[[trait]])
          .coef <- paste0(trait, .levels[-1])
          limma_top <-  as.data.table(topTable(
            fit = limma_fit2,
            coef = .coef,
            number = nrow(beta_matrix),
            adjust.method = "BH",
            p.value = 1,
            sort.by = "none"
          ), keep.rownames = "CpG")[
            j = `:=`(
              "se" = sqrt(limma_fit2[["s2.post"]]) * limma_fit2[["stdev.unscaled"]][, .coef],
              "adj.P.Val" = NULL,
              "B" = NULL,
              "fdr" = p.adjust(P.Value, method = "BH"),
              "contrast" = paste0(trait, ": ", .levels[2], " Vs. ", .levels[1], " (ref)")
            )
          ]
        } else {
          limma_top <-  as.data.table(topTable(
            fit = limma_fit2,
            coef = trait,
            number = nrow(beta_matrix),
            adjust.method = "BH",
            p.value = 1,
            sort.by = "none"
          ), keep.rownames = "CpG")[
            j = `:=`(
              "se" = sqrt(limma_fit2[["s2.post"]]) * limma_fit2[["stdev.unscaled"]][, trait],
              "adj.P.Val" = NULL,
              "B" = NULL,
              "fdr" = p.adjust(P.Value, method = "BH"),
              "contrast" = trait
            )
          ]
        }
      }
      setnames(
        x = limma_top, 
        old = c("logFC", "AveExpr", "t", "P.Value"), 
        new = c("estimate", "avgmvalue_meth", "t_statistic", "pvalue")
      )

      limma_annot <- Reduce(
        f = function(x, y) merge(x, y, by = "CpG", all.x = TRUE), 
        x = list(
          limma_top,
          as.data.table(`names<-`(Locations, paste0("cpg_", names(Locations))), keep.rownames = "CpG"),
          as.data.table(Islands.UCSC, keep.rownames = "CpG"),
          as.data.table(Other, keep.rownames = "CpG")[
            j = .(CpG, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group)
          ]
        )
      )
      
      fwrite(
        x = limma_annot,
        file = file.path(output_directory, trait, glue("{project_name}_EWAS_DMP_{trait}{imodel}.csv.gz"))
      )
      
      # dmr_object <- cpg.annotate(
      #   datatype = "array",
      #   object = beta_matrix[, as.character(pheno[["epic"]])],
      #   what = "M",
      #   arraytype = "EPIC",
      #   analysis.type = "differential",
      #   design = model.matrix(object = form, data = pheno),
      #   coef = 2, 
      #   fdr = 0.05
      # )
      # any_dmr_cpg <- any(dmr_object@ranges@elementMetadata@listData$is.sig)
      # dmr <- as.data.table(suppressPackageStartupMessages({
      #   extractRanges(
      #     dmrcate(
      #       object = dmr_object,
      #       pcutoff = if (any_dmr_cpg) 0.05 else 1,
      #       lambda = 1000, # default
      #       C = NULL, # default
      #       min.cpgs = 5 # default 2 
      #     ),
      #     genome = "hg19"
      #   )
      # }))
      # fwrite(
      #   x = dmr,
      #   file = file.path(
      #     output_directory, trait, 
      #     glue("{project_name}_EWAS_DMR_{trait}{imodel}_{if (any_dmr_cpg) 'fdr005' else 'fdr100'}.csv.gz")
      #   )
      # )
    }
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
