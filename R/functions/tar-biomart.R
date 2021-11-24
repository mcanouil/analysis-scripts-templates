#' get_biomart_information
#' @import data.table
#' @import biomaRt
#' @import httr
get_biomart_information <- function(
  txi,
  ensembl_id,
  rna_level = c("ensembl_gene_id", "ensembl_transcript_id"),
  organism = "hsapiens_gene_ensembl",
  version,
  build = 38
) {
  # "hsapiens_gene_ensembl"
  # "mmusculus_gene_ensembl"
  # "rnorvegicus_gene_ensembl"
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  get_mart <- quote(biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = organism,
    version = version,
    GRCh = if (build == 37) build else NULL
  ))

  mart <- try(eval(get_mart), silent = TRUE)
  if (inherits(mart, "try-error")) mart <- eval(get_mart)
  ensembl_build_version <- sprintf("GRCh%d-%s", build, version)

  if (missing(ensembl_id) || is.null(ensembl_id)) {
    if (any(grepl("counts", names(txi)))) {
      ensembl_id <- rownames(txi[["counts"]])
    } else {
      ensembl_id <- rownames(unlist(txi, recursive = FALSE)[["counts"]])
    }
  }
  list_unique_gene <- list(sub("\\..*$", "", unique(ensembl_id)))

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
  ][j = ensembl_version := ensembl_build_version]

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