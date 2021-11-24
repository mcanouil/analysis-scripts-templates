#' get_biomart_information
#' @import data.table
#' @import biomaRt
get_biomart_information <- function(
  ensembl_id,
  rna_level = c("ensembl_gene_id", "ensembl_transcript_id"),
  organism = "hsapiens_gene_ensembl",
  version
) {
  # "hsapiens_gene_ensembl"
  # "mmusculus_gene_ensembl"
  # "rnorvegicus_gene_ensembl"

  mart <- try(biomaRt::useEnsembl(biomart = "ensembl", dataset = organism, version = version), silent = TRUE)
  if (inherits(mart, "try-error")) {
    mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = organism, version = version)
  }

  list_unique_gene <- list(unique(ensembl_id))

  format_columns <- function(x) {
    out <- paste(setdiff(unique(x), ""), collapse = ";")
    data.table::fifelse(out == "", NA_character_, out)
  }

  ensembl_dt <- data.table::setDT(biomaRt::getBM(
    attributes = c(
      rna_level,
      "chromosome_name",
      "start_position",
      "end_position",
      "external_gene_name"
    ),
    filters = rna_level,
    values = list_unique_gene,
    mart = mart
  ))[
    j = lapply(.SD, format_columns),
    by = rna_level
  ][j = ensembl_version := sprintf("GRCh38-%d", version)]

  entrez_dt <- data.table::setDT(biomaRt::getBM(
    attributes = c(rna_level, "entrezgene_id"),
    filters = rna_level,
    values = list_unique_gene,
    mart = mart
  ))[
    j = lapply(.SD, format_columns),
    by = rna_level
  ]

  uniprot_dt <- data.table::setDT(biomaRt::getBM(
    attributes = c(rna_level, "uniprotswissprot"),
    filters = rna_level,
    values = list_unique_gene,
    mart = mart
  ))[
    j = lapply(.SD, format_columns),
    by = rna_level
  ]

  merge(
    x = ensembl_dt,
    y = merge(x = entrez_dt, y = uniprot_dt, by = rna_level, all = TRUE),
    by = rna_level,
    all.x = TRUE
  )
}