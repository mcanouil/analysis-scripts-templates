### global libraries ===============================================================================
library(targets)
# library(tarchetypes)
# library(here)
# library(data.table)

# targets::tar_renv(extras = "visNetwork", path = "scripts/_dependencies.R")
tar_option_set(cue = tar_cue(mode = "never"))


### project setup ==================================================================================
output_directory <- here("outputs", "ora")
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

threshold_gene_ora <- 0.05
threshold_pathway_ora <- 1


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
  ### ora ------------------------------------------------------------------------------------------
  tar_target(ora,
    command = {
      suppressPackageStartupMessages(library(organism[["go"]], character.only = TRUE))
      list(
        "Reactome" = enrichPathway(
          gene = unique(na.exclude(unlist(strsplit(results[fdr < threshold_gene_ora][["entrezgene_id"]], ";")))),
          universe = unique(na.exclude(unlist(strsplit(results[["entrezgene_id"]], ";")))),
          organism = organism[["reactome"]],
          pvalueCutoff = threshold_pathway_ora,
          pAdjustMethod = "BH",
          readable = TRUE
        ),
        "Gene Ontology Biological Process" = enrichGO(
          gene = unique(results[fdr < threshold_gene_ora][["ensembl_gene_id"]]),
          universe = unique(results[["ensembl_gene_id"]]),
          OrgDb = get(organism[["go"]]),
          keyType = "ENSEMBL",
          ont = "BP",
          pvalueCutoff = threshold_pathway_ora,
          pAdjustMethod = "BH",
          readable = TRUE
        ),
        "Gene Ontology Cellular Component" = enrichGO(
          gene = unique(results[fdr < threshold_gene_ora][["ensembl_gene_id"]]),
          universe = unique(results[["ensembl_gene_id"]]),
          OrgDb = get(organism[["go"]]),
          keyType = "ENSEMBL",
          ont = "CC",
          pvalueCutoff = threshold_pathway_ora,
          pAdjustMethod = "BH",
          readable = TRUE
        ),
        "Gene Ontology Molecular Function" = enrichGO(
          gene = unique(results[fdr < threshold_gene_ora][["ensembl_gene_id"]]),
          universe = unique(results[["ensembl_gene_id"]]),
          OrgDb = get(organism[["go"]]),
          keyType = "ENSEMBL",
          ont = "MF",
          pvalueCutoff = threshold_pathway_ora,
          pAdjustMethod = "BH",
          readable = TRUE
        ),
        "KEGG" = {
          kegg_res <- enrichKEGG(
            gene = unique(na.exclude(unlist(strsplit(results[downstream_pvalue < threshold_gene_ora][["uniprotswissprot"]], ";")))),
            universe = unique(na.exclude(unlist(strsplit(results[["uniprotswissprot"]], ";")))),
            organism = organism[["kegg"]],
            keyType = "uniprot",
            pvalueCutoff = threshold_pathway_ora,
            pAdjustMethod = "BH"
          )
          kegg_res@result <- setDT(kegg_res@result)[
            j = "geneID_symbols" := lapply(
              X = .SD,
              FUN = function(icol) {
                sapply(
                  X = icol, 
                  res = results,
                  FUN = function(x, res) {
                    x <- unlist(setdiff(na.exclude(strsplit(x, "/")), c("", "NA")))
                    if (all(grepl("ENSG", x))) {
                      id <- "ensembl_id"
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
            .SDcols = "geneID"
          ]
          setDF(kegg_res@result)
          kegg_res
        }
      )
    },
    packages = c("data.table", "clusterProfiler", "ReactomePA", organism[["go"]])
  ),
  ### write_ora ------------------------------------------------------------------------------------
  tar_target(write_ora,
    command = {
      write_xlsx(
        x = setNames(lapply(ora, FUN = function(.enrich) {
          if (is.null(.enrich)) return(data.frame())
          .enrich@result
        }), gsub("Gene Ontology", "GO", names(ora))),
        path = file.path(output_directory, "over_representation.xlsx")
      )
      file.path(output_directory, "over_representation.xlsx")
    },
    packages = c("writexl"), 
    format = "file"
  )
)
