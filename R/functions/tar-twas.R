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
      j = `:=`(
        Sample_ID = sub(".genes.results$", "", basename(rnaseq_path))
      )
    ][
      j = `:=`(
        group = factor(sub("-.*", "", Sample_ID), levels = unique(sub("-.*", "", Sample_ID)))
      )
    ]
  } else {

  }
}
#' read_rsem
#' @import tximport
read_rsem <- function(sample_sheet) {
  if (!all(c("rnaseq_path", "Sample_ID") %in% colnames(sample_sheet))) {
    stop("sample_sheet must contain columns \"rnaseq_path\" and \"Sample_ID\"!")
  }
  txi_counts <- tximport::tximport(
    files = setNames(sample_sheet[["rnaseq_path"]], sample_sheet[["Sample_ID"]]),
    type = "rsem",
    txIn = FALSE,
    txOut = FALSE,
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
plot_pca_twas <- function(data, sample_sheet, pca_vars, n_comp = 10, fig_n_comp = 3) {
  n_comp <- min(n_comp, ncol(txi_counts))
  fig_n_comp <- min(fig_n_comp, ncol(txi_counts))

  txi_counts <- data[["counts"]]
  txi_counts <- txi_counts[rowSums(is.na(txi_counts)) == 0, ]

  sample_sheet <- sample_sheet[Sample_ID %in% colnames(txi_counts)]
  if (missing(pca_vars) || is.null(pca_vars)) {
    pca_vars <- colnames(sample_sheet)
  } else {
    pca_vars <- intersect(colnames(sample_sheet), pca_vars)
  }

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

  pca_res <- flashpcaR::flashpca(X = t(data), stand = "sd", ndim = n_comp)

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
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
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
                {
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
                }
            }
          ),
          list(p_inertia)
        ),
        guides = "collect"
      ) +
        patchwork::plot_annotation(
          title = paste0("Structure Detection For: '<i>", ivar, "</i>'"),
          tag_levels = "A",
          theme = ggplot2::theme(plot.title = ggtext::element_markdown())
        )
    })
  )
}
