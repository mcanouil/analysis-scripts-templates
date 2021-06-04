### global libraries ===============================================================================
library(targets)
# library(tarchetypes)
# library(here)
# library(data.table)

# targets::tar_renv(extras = "visNetwork", path = "scripts/_dependencies.R")
tar_option_set(cue = tar_cue(mode = "never"))


### project setup ==================================================================================
output_directory <- here("outputs", "gsea")
dir.create(output_directory, mode = "0775", showWarnings = FALSE)

set_ensembl_version <- 103

organism <- c(
  "reactome" = "human",
  "go" = "org.Hs.eg.db",
  "kegg" = "hsa",
  "ensembl" = "hsapiens_gene_ensembl"
)
# organism <- c(
#   "reactome" = "mouse",
#   "go" = "org.Mm.eg.db",
#   "kegg" = "mmu",
#   "ensembl" = "mmusculus_gene_ensembl"
# )
# organism <- c(
#   "reactome" = "rat",
#   "go" = "org.Rn.eg.db",
#   "kegg" = "rno",
#   "ensembl" = "rnorvegicus_gene_ensembl"
# )

threshold_gene_gsea <- 0.05
threshold_pathway_gsea <- 1


### targets ========================================================================================
list(
  ## results ---------------------------------------------------------------------------------------
  tar_target(results,
    command = {

    },
    packages = c("data.table", "here")
  ),
  ## biomart ---------------------------------------------------------------------------------------
  tar_target(biomart,
    command = {
      set_config(config(ssl_verifypeer = FALSE)) # Fix SSL check error

      mart <- try(useEnsembl(biomart = "ensembl", dataset = organism[["ensembl"]], version = set_ensembl_version), silent = TRUE)
      if (inherits(mart, "try-error")) {
        mart <- useEnsembl(biomart = "ensembl", dataset = organism[["ensembl"]], version = set_ensembl_version)
      }
      
      ensembl_dt <- setDT(getBM(
        attributes = c(
          "ensembl_gene_id",
          "chromosome_name",
          "start_position",
          "end_position",
          "external_gene_name"
        ),
        filters = "ensembl_gene_id",
        values = list(unique(eqtl[["ensembl_gene_id"]])),
        mart = mart
      ))[
        j = lapply(.SD, function(x) {
          out <- paste(setdiff(unique(x), ""), collapse = ";")
          fifelse(out == "", NA_character_, out)
        }),
        by = "ensembl_gene_id"
      ][j = ensembl_version := sprintf("GRCh38-%d", set_ensembl_version)]
      
      entrez_dt <- setDT(getBM(
        attributes = c("ensembl_gene_id", "entrezgene_id"),
        filters = "ensembl_gene_id",
        values = list(unique(eqtl[["ensembl_gene_id"]])),
        mart = mart
      ))[
        j = lapply(.SD, function(x) {
          out <- paste(setdiff(unique(x), ""), collapse = ";")
          fifelse(out == "", NA_character_, out)
        }),
        by = "ensembl_gene_id"
      ]
      
      uniprot_dt <- setDT(getBM(
        attributes = c("ensembl_gene_id", "uniprotswissprot"),
        filters = "ensembl_gene_id",
        values = list(unique(eqtl[["ensembl_gene_id"]])),
        mart = mart
      ))[
        j = lapply(.SD, function(x) {
          out <- paste(setdiff(unique(x), ""), collapse = ";")
          fifelse(out == "", NA_character_, out)
        }),
        by = "ensembl_gene_id"
      ]
      
      merge(
        x = ensembl_dt,
        y = merge(x = entrez_dt, y = uniprot_dt, by = "ensembl_gene_id", all = TRUE),
        by = "ensembl_gene_id",
        all.x = TRUE
      )
    },
    packages = c("biomaRt", "httr", "data.table")
  ),
  ### annot ----------------------------------------------------------------------------------------
  tar_target(annot,
    command = {
      merge(
        x = results,
        y = biomart,
        by = "ensembl_gene_id",
        all.x = TRUE
      )[order(fdr)]
    },
    packages = c("data.table", "here")
  ),
  ### gsea -----------------------------------------------------------------------------------------
  tar_target(gsea,
    command = {
      suppressPackageStartupMessages(library(organism[["go"]], character.only = TRUE))
      list(
        "Reactome" = {
          genes_list <- results[
            downstream_pvalue < threshold_gene_gsea & !is.na(entrezgene_id) & entrezgene_id != ""
          ][
            !duplicated(entrezgene_id)
          ][
            i = order(-slope),
            j = setNames(slope, entrezgene_id)
          ]
          gsePathway(
            geneList = sub(";.*", "", genes_list), 
            organism = organism[["reactome"]],
            pvalueCutoff = threshold_pathway_gsea, 
            pAdjustMethod = "BH"
          )
        },
        "Gene Ontology Biological Process" = {
          genes_list <- results[
            downstream_pvalue < threshold_gene_gsea & !is.na(ensembl_gene_id) & ensembl_gene_id != ""
          ][
            !duplicated(ensembl_gene_id)
          ][
            i = order(-slope),
            j = setNames(slope, ensembl_gene_id)
          ]
          gseGO(
            geneList = sub(";.*", "", genes_list),
            OrgDb = get(organism[["go"]]),
            keyType = "ENSEMBL",
            ont = "BP",
            pvalueCutoff = threshold_pathway_gsea,
            pAdjustMethod = "BH"
          )
        },
        "Gene Ontology Cellular Component" = {
          genes_list <- results[
            downstream_pvalue < threshold_gene_gsea & !is.na(ensembl_gene_id) & ensembl_gene_id != ""
          ][
            !duplicated(ensembl_gene_id)
          ][
            i = order(-slope),
            j = setNames(slope, ensembl_gene_id)
          ]
          gseGO(
            geneList = sub(";.*", "", genes_list),
            OrgDb = get(organism[["go"]]),
            keyType = "ENSEMBL",
            ont = "CC",
            pvalueCutoff = threshold_pathway_gsea,
            pAdjustMethod = "BH"
          )
        },
        "Gene Ontology Molecular Function" = {
          genes_list <- results[
            downstream_pvalue < threshold_gene_gsea & !is.na(ensembl_gene_id) & ensembl_gene_id != ""
          ][
            !duplicated(ensembl_gene_id)
          ][
            i = order(-slope),
            j = setNames(slope, ensembl_gene_id)
          ]
          gseGO(
            geneList = sub(";.*", "", genes_list),
            OrgDb = get(organism[["go"]]),
            keyType = "ENSEMBL",
            ont = "MF",
            pvalueCutoff = threshold_pathway_gsea,
            pAdjustMethod = "BH"
          )
        },
        "KEGG" = {
          genes_list <- results[
            downstream_pvalue < threshold_gene_gsea & !is.na(uniprotswissprot) & uniprotswissprot != ""
          ][
            !duplicated(uniprotswissprot)
          ][
            i = order(-slope),
            j = setNames(slope, uniprotswissprot)
          ]
          gseKEGG(
            geneList = sub(";.*", "", genes_list),
            organism = organism[["kegg"]], 
            keyType = "uniprot", 
            pvalueCutoff = threshold_pathway_gsea, 
            pAdjustMethod = "BH"
          )
        }
      )
    },
    packages = c("data.table", "clusterProfiler", "ReactomePA", organism[["go"]])
  ),
  ### write_gsea -----------------------------------------------------------------------------------
  tar_target(write_gsea,
    command = {
      write_xlsx(
        x = setNames(lapply(gsea, FUN = function(.enrich) {
          if (is.null(.enrich) || nrow(.enrich@result) == 0) return(data.frame())
          merge(
            x = setDT(.enrich@result), 
            y = setnames(as.data.table(
              x = sapply(.enrich@geneSets, function(.l) {
                paste(sort(intersect(.l, names(.enrich@geneList))), collapse = "/")
              }), 
              keep.rownames = TRUE
            ), c("ID", "genes_set")),
            by = "ID"
          )[
            j = peripheral_enrichment := {
              peripheral_set <- setdiff(
                unlist(tstrsplit(genes_set, "/"), recursive = TRUE),
                unlist(tstrsplit(core_enrichment, "/"), recursive = TRUE)
              )
              if (length(peripheral_set) == 0) {
                NA_character_
              } else  {
                paste(peripheral_set, collapse = "/")
              }
            },
            by = "ID"
          ][
            j = c("core_enrichment_symbols", "genes_set_symbols", "peripheral_enrichment_symbols") := lapply(
              X = .SD,
              FUN = function(icol) {
                sapply(
                  X = icol, 
                  res = results,
                  FUN = function(x, res) {
                    x <- unlist(setdiff(na.exclude(strsplit(x, "/")), c("", "NA")))
                    if (all(grepl("ENSG", x))) {
                      id <- "ensembl_gene_id"
                    } else if (
                      all(grepl("[[:alpha:]]", substr(x, 1, 1)) & 
                        !grepl("[[:digit:]]", substr(x, 1, 1)))
                    ) {
                      id <- "uniprotswissprot"
                    } else {
                      id <- "entrezgene_id"
                    }
                    gene_symbols <- unname(setNames(res[["external_gene_name"]], res[[id]])[x])
                    
                    if (length(gene_symbols) == 0) return(NA_character_)
                    
                    paste(gene_symbols, collapse = "/")
                  }
                )
              }
            ),
            .SDcols = c("core_enrichment", "genes_set", "peripheral_enrichment")
          ]
        }), gsub("Gene Ontology", "GO", names(gsea))), 
        path = file.path(output_directory, "gene_set_enrichment.xlsx")
      )
      file.path(output_directory, "gene_set_enrichment.xlsx")
    },
    packages = c("writexl"), 
    format = "file"
  )
)
