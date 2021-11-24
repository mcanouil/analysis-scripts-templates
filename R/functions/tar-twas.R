#' qc_sample_sheet_twas
#' @import data.table
qc_sample_sheet_twas <- function(phenotype, run_path) {
  if (missing(phenotype) || is.null(phenotype)) {
    data.table::data.table(
      rnaseq_path = list.files(
        path = file.path(run_path, "Output", "RSEM"),
        pattern = ".genes.results$",
        full.names = TRUE
      )
    )[
      j = `:=`(Sample_ID = sub(".genes.results$", "", basename(rnaseq_path)))
    ][
      j = `:=`(
        group = factor(x = sub("-.*", "", Sample_ID), levels = c("bsa", "palmitate"))
      )
    ][
      j = `:=`(rep = sub("^[^-]+-(r[^-]+).*", "\\1", Sample_ID))
    ][
      j = `:=`(group.rep = factor(sprintf("%s.%s", group, rep)))
    ]
  } else {

  }
}

#' read_rsem
#' @import tximport
read_rsem <- function(sample_sheet, rna_level = c("ensembl_gene_id", "ensembl_transcript_id")) {
  rna_level <- c("ensembl_gene_id" = "genes", "ensembl_transcript_id" = "isoforms")[rna_level][1]
  if (!all(c("rnaseq_path", "Sample_ID") %in% colnames(sample_sheet))) {
    stop("sample_sheet must contain columns \"rnaseq_path\" and \"Sample_ID\"!")
  }
  rsem_files <- `names<-`(sample_sheet[["rnaseq_path"]], sample_sheet[["Sample_ID"]])

  txi_counts <- tximport::tximport(
    files = sub("\\.genes\\.results$", sprintf(".%s.results", rna_level), rsem_files),
    type = "rsem",
    txIn = rna_level == "transcript",
    txOut = rna_level == "transcript",
    countsFromAbundance = "no"
  )
  txi_counts$length[txi_counts$length == 0] <- 1

  txi_counts
}

#' plot_pca_twas
#' @import data.table
#' @import flashpcaR
#' @import ggplot2
#' @import ggtext
#' @import patchwork
#' @import scales
#' @import stats
#' @import utils
#' @import DESeq2
#' @import MatrixGenerics
plot_pca_twas <- function(txi, sample_sheet, pca_vars, n_comp = 10, fig_n_comp = 3) {
  if (missing(pca_vars) || is.null(pca_vars)) {
    pca_vars <- colnames(sample_sheet)
  } else {
    pca_vars <- intersect(colnames(sample_sheet), pca_vars)
  }

  dds_pheno <- na.exclude(
    sample_sheet[
      j = .SD,
      .SDcols = unique(c("Sample_ID", pca_vars))
    ]
  )

  dds_counts <- lapply(
    X = txi,
    .sample = as.character(dds_pheno[["Sample_ID"]]),
    FUN = function(.l, .samples) if (is.matrix(.l)) .l[, .samples] else .l
  )

  dds <- DESeq2::DESeqDataSetFromTximport(txi = dds_counts, colData = dds_pheno, design = ~ 1)
  dds <- dds[
    MatrixGenerics::rowVars(DESeq2::counts(dds)) != 0 &
      rowMeans(DESeq2::counts(dds)) > 1 &
      MatrixGenerics::rowMedians(DESeq2::counts(dds)) > 0,
  ]

  txi_counts <- DESeq2::counts(dds)
  sample_sheet <- sample_sheet[Sample_ID %in% colnames(txi_counts)]

  keep_technical <- names(which(sapply(sample_sheet[
    j = lapply(.SD, function(x) {
      (data.table::uniqueN(x) > 1 & data.table::uniqueN(x) < length(x)) | is.numeric(x)
    }),
    .SDcols = pca_vars
  ], isTRUE)))

  variables_excluded <- setdiff(pca_vars, keep_technical)
  if (length(variables_excluded) != 0) {
    message(paste(
      "The following variables have been excluded (null variances or confounding with samples):\n",
      paste("+", variables_excluded),
      "\n",
      sep = "\n"
    ))
  }
  if (length(keep_technical) == 0) return(NULL)

  n_comp <- min(n_comp, ncol(txi_counts) - 5)
  fig_n_comp <- min(fig_n_comp, ncol(txi_counts) - 5)
  pca_res <- flashpcaR::flashpca(X = t(txi_counts), stand = "sd", ndim = n_comp)

  pca_dfxy <- data.table::as.data.table(pca_res[["vectors"]], keep.rownames = "Sample_ID")
  data.table::setnames(
    x = pca_dfxy,
    old = setdiff(names(pca_dfxy), "Sample_ID"),
    new = sprintf("PC%02d", as.numeric(gsub("V", "", setdiff(names(pca_dfxy), "Sample_ID"))))
  )
  pca_dfxy <- merge(x = pca_dfxy, y = sample_sheet, by = "Sample_ID")

  p_inertia <- ggplot2::ggplot(
    data = data.table::data.table(
      y = pca_res[["pve"]],
      x = sprintf("PC%02d", seq_along(pca_res[["pve"]]))
    )[x %in% sprintf("PC%02d", seq_len(fig_n_comp))]
  ) +
    ggplot2::aes(
      x = paste0(
        x,
        "<br><i style='font-size:8pt;'>(",
        scales::percent_format(accuracy = 0.01, suffix = " %")(y),
        ")</i>"
      ),
      y = y
    ) +
    ggplot2::geom_col(
      width = 1,
      colour = "white",
      fill = scales::viridis_pal(begin = 0.5, end = 0.5)(1),
      na.rm = TRUE
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1, suffix = " %"),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      x = "Principal Components",
      y = "Variance Explained"
    ) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggtext::element_markdown(),
      plot.subtitle = ggtext::element_markdown(face = "italic"),
      axis.text.x = ggtext::element_markdown()
    )

  asso_dt <- data.table::melt(
    data = pca_dfxy,
    measure.vars = grep("^PC[0-9]+$", names(pca_dfxy), value = TRUE),
    variable.name = "pc",
    value.name = "values"
  )[pc %in% sprintf("PC%02d", seq_len(n_comp))][
    j = {
      m <- stats::model.matrix(
        object = stats::as.formula(
          object = paste0("values ~ ", paste(sprintf("`%s`", keep_technical), collapse = " + "))
        ),
        data = .SD
      )

      if (qr(m)$rank == ncol(m)) {
        out <- data.table::as.data.table(
          stats::anova(
            stats::lm(
              formula = stats::as.formula(
                object = paste0("values ~ ", paste(sprintf("`%s`", keep_technical), collapse = " + "))
              ),
              data = .SD
            )
          ),
          keep.rownames = "term"
        )[term != "Residuals"]
      } else {
        out <- data.table::rbindlist(
          lapply(X = sprintf("`%s`", keep_technical), .data = .SD, FUN = function(.x, .data) {
            data.table::as.data.table(
              stats::anova(
                stats::lm(
                  formula = stats::as.formula(paste0("values ~ ", .x)),
                  data = .SD
                )
              ),
              keep.rownames = "term"
            )[term != "Residuals"]
          })
        )
      }
      out[j = full_rank := qr(m)$rank == ncol(m)][j = term := gsub("`", "", term)]
    },
    by = "pc"
  ]

  p_association <- ggplot2::ggplot(data = asso_dt) +
    ggplot2::aes(
      x = factor(.data[["pc"]]),
      y = factor(
        x = .data[["term"]],
        levels = data.table::setorderv(
          x = data.table::dcast(
            data = asso_dt[j = list(pc, term, `Pr(>F)` = data.table::fifelse(`Pr(>F)` <= 0.1, `Pr(>F)`, NA_real_))],
            formula = term ~ pc,
            value.var = "Pr(>F)"
          ),
          cols = levels(asso_dt[["pc"]])[seq_len(n_comp)],
          order = -1
        )[["term"]]
      ),
      fill = .data[["Pr(>F)"]]
    ) +
    ggplot2::geom_tile(colour = "white", na.rm = TRUE) +
    ggtext::geom_richtext(
      mapping = ggplot2::aes(
        label = gsub(
          pattern = "(.*)e([-+]*)0*(.*)",
          replacement = "\\1<br>&times;<br>10<sup>\\2\\3</sup>",
          x = format(.data[["Pr(>F)"]], digits = 2, nsmall = 2, scientific = TRUE)
        )
      ),
      colour = "white",
      fill = NA,
      label.colour = NA,
      size = 4,
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_viridis_c(na.value = "grey85", end = 0.75, limits = c(0, 0.1)) +
    ggplot2::scale_x_discrete(
      expand = c(0, 0),
      labels = function(x) {
        paste0(
          x,
          "<br><i style='font-size:8pt;'>(",
          format(
            x = pca_res[["pve"]][as.numeric(gsub("PC", "", x))] * 100,
            digits = 2,
            nsmall = 2
          ),
          " %)</i>"
        )
      }
    ) +
    ggplot2::scale_y_discrete(expand = c(0, 0), labels = toupper) +
    ggplot2::labs(
      x = "Principal Components",
      y = "Variables",
      title = "Association Tests Between Variables And Principal Components",
      caption = ifelse(
        test = all(asso_dt[["full_rank"]]),
        yes = "Variables are tested against principal components using ANOVA.",
        no = paste(
          "Variables are independently tested against principal components using ANOVA",
          "(*i.e.*, model matrix is not full rank)."
        )
      ),
      fill = "P-Value"
    ) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggtext::element_markdown(),
      plot.subtitle = ggtext::element_markdown(face = "italic"),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggtext::element_markdown()
    )

  c(
    p_association = list(p_association),
    lapply(stats::setNames(keep_technical, keep_technical), function(ivar) {
      patchwork::wrap_plots(
        c(
          apply(
            X = utils::combn(sprintf("PC%02d", seq_len(fig_n_comp)), 2),
            MARGIN = 2,
            FUN = function(x) {
              ggplot2::ggplot(data = pca_dfxy[j = .SD, .SDcols = c(ivar, x)]) +
                ggplot2::aes(x = .data[[x[1]]], y = .data[[x[2]]], colour = .data[[ivar]]) +
                ggplot2::geom_hline(yintercept = 0, linetype = 2, na.rm = TRUE) +
                ggplot2::geom_vline(xintercept = 0, linetype = 2, na.rm = TRUE) +
                ggplot2::geom_point(na.rm = TRUE) +
                (
                  if (is.numeric(pca_dfxy[[ivar]])) {
                    ggplot2::scale_colour_viridis_c(
                      name = NULL,
                      begin = 0,
                      end = 0.75
                    )
                  } else {
                    list(
                      ggplot2::stat_ellipse(type = "norm", na.rm = TRUE, show.legend = FALSE),
                      ggplot2::scale_colour_viridis_d(
                        name = NULL,
                        begin = if (pca_dfxy[j = data.table::uniqueN(.SD), .SDcols = ivar] == 2) 0.25 else 0,
                        end = 0.75,
                        guide = ggplot2::guide_legend(override.aes = list(size = 4))
                      ),
                      if (length(unique(pca_dfxy[[ivar]])) > 10) {
                        ggplot2::theme(legend.position = "none")
                      } else {
                        NULL
                      }
                    )
                  }
                ) +
                ggplot2::theme(
                  plot.title.position = "plot",
                  plot.caption.position = "plot",
                  plot.title = ggtext::element_markdown(),
                  plot.subtitle = ggtext::element_markdown(face = "italic"),
                  axis.text.x = ggtext::element_markdown()
                )
            }
          ),
          list(p_inertia)
        ),
        guides = "collect"
      ) +
        patchwork::plot_annotation(
          title = sprintf("Structure Detection For: '\"<i>%s</i>\"", ivar),
          tag_levels = "A",
          theme = ggplot2::theme(plot.title = ggtext::element_markdown())
        )
    })
  )
}

#' do_twas
#' @import data.table
#' @import DESeq2
#' @import S4Vectors
#' @import MatrixGenerics
#' @import utils
#' @import stats
do_twas <- function(txi, sample_sheet, model, path, rna_level = c("ensembl_gene_id", "ensembl_transcript_id"), biomart) {
  rna_level <- rna_level[1]
  if (is.null(model[["covariates"]]) || nchar(model[["covariates"]]) == 0) {
    covariates <- NULL
  } else {
    covariates <- all.vars(stats::as.formula(sprintf(" ~ %s", model[["covariates"]])))
  }

  form <- stats::as.formula(paste0("~ ", paste(c(model[["raw_trait"]], covariates), collapse = " + ")))

  pheno_dt <- sample_sheet[
    j = na.exclude(.SD),
    .SDcols = c("Sample_ID", all.vars(form))
  ]

  raw_trait <- all.vars(stats::as.formula(paste0("~", model[["raw_trait"]])))
  if (grepl("factor\\(.+\\)", model[["raw_trait"]])) {
    pheno_dt[
      j = c(raw_trait) := lapply(
        X = .SD,
        FUN = function(x) {
          eval(parse(text = sub(raw_trait, "x", model[["raw_trait"]])))
        }
      ),
      .SDcols = c(raw_trait)
    ]
    model[["raw_trait"]] <- raw_trait
    form <- stats::as.formula(paste0("~ ", paste(c(model[["raw_trait"]], covariates), collapse = " + ")))
  }

  trait_values <- pheno_dt[[model[["raw_trait"]]]]

  dds_counts <- lapply(
    X = txi,
    .sample = as.character(pheno_dt[["Sample_ID"]]),
    FUN = function(.l, .samples) if (is.matrix(.l)) .l[, .samples] else .l
  )

  message("Performing DESeq2 regression ...")

  dds <- DESeq2::DESeqDataSetFromTximport(txi = dds_counts, colData = pheno_dt, design = form)
  dds <- dds[
    MatrixGenerics::rowVars(DESeq2::counts(dds)) != 0 &
      rowMeans(DESeq2::counts(dds)) > 1 &
      MatrixGenerics::rowMedians(DESeq2::counts(dds)) > 0,
  ]
  stats_dds <- DESeq2::replaceOutliers(
    object = DESeq2::nbinomWaldTest(
      object = DESeq2::estimateDispersions(DESeq2::estimateSizeFactors(dds)),
      maxit = 1000
    ),
    minReplicates = ncol(dds)
  )

  is_issue <- data.table::data.table(
    x = rownames(dds),
    converge = S4Vectors::mcols(stats_dds)$betaConv
  )
  data.table::setnames(x = is_issue, old = "x", new = rna_level)

  if (is.factor(trait_values)) {
    results_dt <- data.table::rbindlist(apply(
      X = utils::combn(levels(trait_values), 2),
      MARGIN = 2,
      .trait = raw_trait,
      .dds = stats_dds,
      FUN = function(.contrast, .trait, .dds) {
        results_dds <- DESeq2::results(
          object = .dds,
          contrast = c(.trait, .contrast[2], .contrast[1]),
          pAdjustMethod = "BH",
          independentFiltering = FALSE,
          cooksCutoff = FALSE
        )
        results_dt <- data.table::as.data.table(results_dds, keep.rownames = rna_level)
        data.table::setnames(results_dt, old = "padj", new = "fdr")
        results_dt[
          j = `:=`("contrast" = sprintf("%s: %s Vs. %s (ref)", .trait, .contrast[2], .contrast[1]))
        ]
        results_dt
      }
    ))
  } else {
    results_dds <- DESeq2::results(
      object = stats_dds,
      name = raw_trait,
      pAdjustMethod = "BH",
      independentFiltering = FALSE,
      cooksCutoff = FALSE
    )
    results_dt <- data.table::as.data.table(results_dds, keep.rownames = rna_level)
    data.table::setnames(results_dt, old = "padj", new = "fdr")
    results_dt[j = `:=`("contrast" = raw_trait)]
  }

  results_dt[j = `:=`("Trait" = raw_trait, "n" = ncol(stats_dds))]

  results_avg_tpm <- Reduce(
    f = function(x, y) merge(x, y, by = rna_level, all.x = TRUE),
    x = lapply(
      X = levels(trait_values),
      FUN = function(gl) {
        data.table::setnames(
          x = data.table::as.data.table(
            x = rowMeans(dds_counts[["abundance"]][, trait_values %in% gl]),
            keep.rownames = TRUE
          ),
          new = c(rna_level, sprintf("TPM_%s", gl))
        )
      }
    ),
    init = merge(x = results_dt, y = is_issue, by = rna_level)
  )

  results_file <- sprintf("%s/twas_%s_%s.csv.gz", path, model[["raw_trait"]], model[["tar_group"]])
  dir.create(
    path = dirname(results_file),
    recursive = TRUE,
    mode = "0755",
    showWarnings = FALSE
  )

  if (is.null(biomart) || missing(biomart)) {
    data.table::fwrite(x = results_avg_tpm[order(fdr)], file = results_file)
  } else {
    data.table::fwrite(
      x =  merge(x = results_avg_tpm, y = biomart, by = rna_level, all.x = TRUE)[order(fdr)],
      file = results_file
    )
  }

  message(sprintf("Writing results to \"%s\"!", results_file))

  results_file
}