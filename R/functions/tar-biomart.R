#' get_biomart_information
#' @import data.table
#' @import biomaRt
#' @import httr
get_biomart_information <- function(
  ensembl_id,
  rna_level = c("ensembl_gene_id", "ensembl_transcript_id"),
  organism = c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl", "rnorvegicus_gene_ensembl"),
  version,
  build = 38
) {
  rna_level <- match.arg(rna_level, c("ensembl_gene_id", "ensembl_transcript_id"))
  organism <- match.arg(organism, c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl", "rnorvegicus_gene_ensembl"))
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

  if (!is.null(biomaRt::searchAttributes(mart, pattern = "entrezgene_id"))) {
    entrez_dt <- data.table::setDT(biomaRt::getBM(
      attributes = c(rna_level, "entrezgene_id"),
      filters = rna_level,
      values = list_unique_gene,
      mart = mart
    ))[
      j = lapply(.SD, format_columns),
      by = rna_level
    ]
  }

  if (!is.null(biomaRt::searchAttributes(mart, pattern = "uniprotswissprot"))) {
    uniprot_dt <- data.table::setDT(biomaRt::getBM(
      attributes = c(rna_level, "uniprotswissprot"),
      filters = rna_level,
      values = list_unique_gene,
      mart = mart
    ))[
      j = lapply(.SD, format_columns),
      by = rna_level
    ]
  }

  datasets_exists <- sapply(c("entrez_dt", "uniprot_dt"), exists)

  if (all(datasets_exists)) {
    merge(
      x = ensembl_dt,
      y = merge(x = entrez_dt, y = uniprot_dt, by = rna_level, all = TRUE),
      by = rna_level,
      all.x = TRUE
    )
  } else {
    merge(
      x = ensembl_dt,
      y = get(names(which(datasets_exists))),
      by = rna_level,
      all.x = TRUE
    )
  }
}
