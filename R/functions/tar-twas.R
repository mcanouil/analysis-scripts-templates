#' qc_sample_sheet_twas
#' @import data.table
qc_sample_sheet_twas <- function(phenotype, run_path) {
  if (missing(phenotype) || is.null(phenotype))
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

#' do_pca
#' @import data.table
#' @import flashpcaR
#' @import patchwork
#' @import scales
#' @import ggplot2
#' @import ggtext
#' @import ragg
do_pca <- function() {
  n_comp <- 5
  fig_n_comp <- 3
  keep_technical <- c(
    "SBA Duration (Weeks)" = "weeks", 
    "Group" = "group", 
    "Quantity of RNA (ng)" = "Quantity of RNA (ng)", 
    "Amount (µl)" = "Amount (µl)", 
    "260/280" = "260/280"
  )
  
  dir.create(file.path(output_directory, "pca"), showWarnings = FALSE)
  
  pca_data <- txi[["counts"]]
  pca_phenotypes <- sample_sheet[Sample_ID %in% colnames(pca_data)]
  pca_res <- flashpca(X = t(pca_data), stand = "sd", ndim = n_comp)
  pca_dfxy <- as.data.table(pca_res[["vectors"]], keep.rownames = "Sample_ID")
  setnames(
    x = pca_dfxy, 
    old = setdiff(names(pca_dfxy), "Sample_ID"), 
    new = sprintf("PC%02d", as.numeric(gsub("V", "", setdiff(names(pca_dfxy), "Sample_ID"))))
  )
  pca_dfxy <- merge(x = pca_dfxy, y = pca_phenotypes, by = "Sample_ID")
  p_inertia <- ggplot(
    data = data.table(
      y = pca_res[["pve"]],
      x = sprintf("PC%02d", seq_along(pca_res[["pve"]]))
    )[x %in% sprintf("PC%02d", 1:fig_n_comp)][
      j = x := paste0(x, "<br><i style='font-size:6pt;'>(", percent_format(accuracy = 0.01, suffix = " %")(y), ")</i>")
    ]
  ) +
    aes(x = .data[["x"]], y = .data[["y"]]) +
    geom_col(
      width = 1,
      colour = "white",
      fill = viridis_pal(begin = 0.5, end = 0.5)(1),
      na.rm = TRUE
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 0.1, suffix = " %"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = "Principal Components",
      y = "Contribution"
    ) +
    theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
      plot.caption = element_markdown(face = "italic", size = rel(0.65)),
      axis.title.x = element_markdown(),
      axis.text.x = element_markdown(),
      axis.title.y = element_markdown(),
      axis.text.y = element_markdown(),
      panel.grid.minor = element_blank()
    ) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

  asso_dt <- melt(
    data = pca_dfxy,
    measure.vars = grep("^PC[0-9]+$", names(pca_dfxy), value = TRUE),
    variable.name = "pc",
    value.name = "values"
  )[pc %in% sprintf("PC%02d", 1:n_comp)][,
    {
      m <- model.matrix(
        object = as.formula(
          object = paste0("values ~ ", paste(paste0("`", keep_technical, "`"), collapse = " + "))
        ),
        data = .SD
      )

      if (qr(m)$rank == ncol(m)) {
        out <- as.data.table(
          anova(
            lm(
              formula = as.formula(
                object = paste0("values ~ ", paste(paste0("`", keep_technical, "`"), collapse = " + "))
              ),
              data = .SD
            )
          ),
          keep.rownames = "term"
        )[term != "Residuals"]
      } else {
        out <- rbindlist(
          lapply(X = keep_technical, .data = .SD, FUN = function(.x, .data) {
            as.data.table(
              anova(
                lm(
                  formula = as.formula(paste0("values ~ `", .x, "`")),
                  data = .SD
                )
              ),
              keep.rownames = "term"
            )[term != "Residuals"]
          })
        )
      }
      out[, full_rank := qr(m)$rank == ncol(m)]
    },
    by = "pc"
  ]

  pca_asso <- ggplot(data = asso_dt) +
    aes(
      x = factor(.data[["pc"]]),
      y = factor(
        x = .data[["term"]],
        levels = setorderv(
          x = dcast(
            data = asso_dt[j = list(pc, term, `Pr(>F)` = fifelse(`Pr(>F)` <= 0.1, `Pr(>F)`, NA_real_))],
            formula = term ~ pc,
            value.var = "Pr(>F)"
          ),
          cols = levels(asso_dt[["pc"]])[1:n_comp],
          order = -1
        )[["term"]]
      ),
      fill = .data[["Pr(>F)"]]
    ) +
    geom_tile(colour = "white", na.rm = TRUE) +
    geom_richtext(
      mapping = aes(
        label = (function(x) {
          x_fmt <- x
          x_fmt[as.numeric(x) < 0.01] <- gsub(
            pattern = "(.*)e([-+]*)0*(.*)",
            replacement = "\\1<br>&times;<br>10<sup>\\2\\3</sup>",
            x = format(x_fmt[as.numeric(x) < 0.01], digits = 2, nsmall = 2, scientific = TRUE)
          )
          x_fmt[as.numeric(x) >= 0.01] <- format(
            x = as.numeric(x_fmt[as.numeric(x) >= 0.01]), 
            digits = 2, nsmall = 2
          )
          x_fmt
        })(.data[["Pr(>F)"]])
      ),
      colour = "white",
      fill = NA,
      label.colour = NA,
      size = 2.5,
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(end = 0.75, limits = c(0, 0.1)) +
    theme(panel.grid = element_blank()) +
    scale_x_discrete(
      expand = c(0, 0),
      labels = function(x) {
        paste0(
          x, "<br><i style='font-size:5pt;'>(",
          format(
            x = pca_res[["pve"]][as.numeric(gsub("PC", "", x))] * 100,
            digits = 2,
            nsmall = 3
          ),
          " %)</i>"
        )
      }
    ) +
    scale_y_discrete(
      expand = c(0, 0),
      labels = function(x) names(keep_technical[match(gsub("`", "", x), keep_technical)])
    ) +
    labs(
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
    theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
      plot.caption = element_markdown(face = "italic", size = rel(0.65)),
      axis.title.x = element_markdown(),
      axis.text.x = element_markdown(),
      axis.title.y = element_markdown(),
      axis.text.y = element_markdown(),
      panel.grid.minor = element_blank()
    ) +
    theme(plot.caption = element_markdown())

  agg_png(
    filename = file.path(output_directory, "pca", "gene_pca_asso.png"),
    width = 16, height = 12, units = "cm", res = 300, scaling = 1
  )
    print(pca_asso)
  invisible(dev.off())

  for (ivar in keep_technical) {
    p <- wrap_plots(
      c(
        plot_planes(pca_dfxy, ivar, fig_n_comp),
        list(p_inertia)
      ),
      guides = "collect"
    ) +
      plot_annotation(
        title = paste0(
          "Structure Detection For: '<i>",
          names(keep_technical[match(ivar, keep_technical)]),
          "</i>'"
        ),
        tag_levels = "A",
        theme = theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
        theme(
          plot.title.position = "plot",
          plot.caption.position = "plot",
          plot.title = element_markdown(),
          plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
          plot.caption = element_markdown(face = "italic", size = rel(0.65)),
          axis.title.x = element_markdown(),
          axis.text.x = element_markdown(),
          axis.title.y = element_markdown(),
          axis.text.y = element_markdown(),
          panel.grid.minor = element_blank()
        )
      )

    agg_png(
      filename = sprintf("%s/pca/gene_pca_planes_%s.png",
        output_directory, 
        tolower(gsub("[.]+$|^[.]+|^X+", "", gsub("[.]+", ".", make.names(ivar))))
      ),
      width = 16, height = 12, units = "cm", res = 300, scaling = 1
    )
      print(p)
    invisible(dev.off())
  }

  c(
    file.path(output_directory, "pca", "gene_pca_asso.png"),
    sprintf("%s/pca/gene_pca_planes_%s.png",
      output_directory, 
      tolower(gsub("[.]+$|^[.]+|^X+", "", gsub("[.]+", ".", make.names(ivar))))
    )
  )
}

#' do_twas
#' @import data.table
#' @import DESeq2
do_twas <- function() {
  fmodel <- as.formula(sprintf("~ %s", model))
  dds_pheno <- na.exclude(
    sample_sheet[
      i = weeks %in% eval(str2expression(samples)), 
      j = .SD, 
      .SDcols = c("Sample_ID", all.vars(fmodel))
    ]
  )
  dds_counts <- lapply(
    X = txi,
    .sample = as.character(dds_pheno[["Sample_ID"]]),
    FUN = function(.l, .samples) if (is.matrix(.l)) .l[, .samples] else .l
  )

  dds <- DESeqDataSetFromTximport(txi = dds_counts, colData = dds_pheno, design = fmodel)
  dds <- dds[rowVars(counts(dds)) != 0 & rowMeans(counts(dds)) > 1 & rowMedians(counts(dds)) > 0, ]
  stats_dds <- replaceOutliers(
    object = nbinomWaldTest(estimateDispersions(estimateSizeFactors(dds)), maxit = 1000),
    minReplicates = 3 # ncol(dds)
  )

  is_issue <- data.table(
    x = rownames(dds),
    converge = mcols(stats_dds)$betaConv
  )
  setnames(x = is_issue, old = "x", new = "ensembl_gene_id")

  results_dt <- rbindlist(lapply(
    X = resultsNames(stats_dds)[-1],
    .data = stats_dds,
    FUN = function(.trait, .data) {
      setnames(
        x = as.data.table(
          x = results(
            object = .data,
            name = .trait,
            pAdjustMethod = "BH",
            independentFiltering = TRUE,
            cooksCutoff = TRUE
          ),
          keep.rownames = "ensembl_gene_id"
        )[j = contrast := .trait],
        old = "padj",
        new = "fdr"
      )
    }
  ))

  results_annot_dt <- merge(
    x = merge(x = results_dt, y = is_issue, by = "ensembl_gene_id"),
    y = biomart,
    by = "ensembl_gene_id",
    all.x = TRUE
  )
  
  Reduce(
    f = function(x, y) merge(x, y, by = "ensembl_gene_id", all.x = TRUE),
    x = lapply(
      X = levels(dds_pheno[["group"]]),
      FUN = function(gl) {
        setnames(
          x = as.data.table(
            x = rowMeans(dds_counts[["abundance"]][, dds_pheno[["group"]] %in% gl]), 
            keep.rownames = TRUE
          ),
          new = c("ensembl_gene_id", sprintf("TPM_%s", gl))
        )
      }
    ), 
    init = results_annot_dt
  )
}

#' plot_planes_twas
#' @import ggplot2
#' @import ggtext
#' @import data.table
#' @import utils
plot_planes_twas <- function(pca_dfxy, ivar, fig_n_comp) {
  apply(
    X = utils::combn(sprintf("PC%02d", seq_len(fig_n_comp)), 2),
    MARGIN = 2,
    FUN = function(.x) {
      ggplot2::ggplot(data = pca_dfxy[, .SD, .SDcols = c(ivar, .x)]) +
        ggplot2::aes(x = .data[[.x[1]]], y = .data[[.x[2]]], colour = .data[[ivar]]) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, na.rm = TRUE) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2, na.rm = TRUE) +
        ggplot2::geom_point(size = 0.75, na.rm = TRUE) +
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
                begin = if (pca_dfxy[, data.table::uniqueN(.SD), .SDcols = ivar] == 2) 0.25 else 0,
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
        } +
        ggplot2::theme(
          plot.title.position = "plot",
          plot.caption.position = "plot",
          plot.title = ggtext::element_markdown(),
          plot.subtitle = ggtext::element_markdown(face = "italic"),
          axis.text.y = ggtext::element_markdown(),
          legend.position = "none"
        )
    }
  )
}


#' plot_volcano_twas
#' @import ggplot2
#' @import ggtext
#' @import data.table
#' @import ggrepel
#' @import stats
plot_volcano_twas <- function(file, model) {
  raw_trait <- all.vars(stats::as.formula(paste0("~", model[["raw_trait"]])))

  dt <- data.table::fread(file)[
    j = gene_label_min := data.table::fifelse(
      test = pvalue == min(pvalue, na.rm = TRUE) &
        !is.na(UCSC_RefGene_Name) &
        !UCSC_RefGene_Name %in% c("", "NA"),
      yes = paste(unique(unlist(strsplit(gsub(",", ";", UCSC_RefGene_Name), ";"))), collapse = ";"),
      no = NA_character_
    ),
    by = "UCSC_RefGene_Name"
  ][
    i = pvalue > 0.05,
    j = pvalue := NA_real_
  ][
    j = file := basename(file)
  ][
    j = cpg_chr := sub("chr", "", cpg_chr)
  ][
    j = c("estimate", "cpg_chr", "cpg_pos", "pvalue", "gene_label_min")
  ][order(pvalue)]

  if (is.numeric(dt[["cpg_chr"]])) {
    dt[j = "cpg_chr" := lapply(.SD, as.character), .SDcols = "cpg_chr"]
  }

  if (dt[!is.na(gene_label_min), .N] > 10) {
    dt[which(!is.na(gene_label_min))[-c(1:10)], gene_label_min := NA_character_]
  }

  alpha <- 0.05 / nrow(dt)

  ggplot2::ggplot(dt) +
    ggplot2::aes(
      x = .data[["estimate"]],
      y = .data[["pvalue"]],
      colour = abs(.data[["estimate"]])
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point(size = 0.60) +
    ggplot2::annotate(
      geom = "rect",
      xmin = -Inf, xmax = Inf, ymin = 1, ymax = alpha,
      fill = "#b22222", alpha = 0.2, colour = NA,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = alpha, linetype = 2, colour = "#b22222") +
    ggplot2::scale_colour_viridis_c(trans = "sqrt", limits = c(0, NA)) +
    ggplot2::scale_y_continuous(
      trans = pval_trans(md = TRUE),
      expand = ggplot2::expansion(mult = c(0, 0.2))
    ) +
    ggplot2::coord_cartesian(ylim = c(0.05, NA)) +
    ggrepel::geom_label_repel(
      mapping = ggplot2::aes(label = .data[["gene_label_min"]]),
      show.legend = FALSE,
      min.segment.length = 0,
      size = 2.5,
      na.rm = TRUE,
      max.overlaps = Inf
    ) +
    ggplot2::labs(
      x = "Estimate",
      y = "P-value",
      colour = "Estimate",
      title = model[["pretty_trait"]],
      subtitle = toupper(paste(raw_trait, "= ", model[["covariates"]]))
    ) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggtext::element_markdown(),
      plot.subtitle = ggtext::element_markdown(face = "italic"),
      axis.text.y = ggtext::element_markdown(),
      legend.position = "none"
    )
}

# draw_volcano <- function(
#   data, x, y,
#   label_x = "Fold-Change (log<sub>2</sub>)",
#   label_y = "P-value",
#   label_colour = label_x,
#   alpha = 0.05,
#   max_p = 1
# ) {
#   data <- data.table::as.data.table(data)#[, .SD, .SDcols = c(x, y)]
#   data.table::setnames(data, y, "pvalue", skip_absent = TRUE)
#   if (!"contrast" %in% names(data)) data[j = contrast := "Volcano plot"]

#   if ("OR" %in% x) {
#     label_x <- gsub("Fold Change", "OR", label_x)
#     data[j = c(x) := log2(.SD), .SDcols = x]
#   }

#   out <- lapply(unique(data[["contrast"]]), function(.contrast) {
#     ggplot2::ggplot(data[contrast %in% .contrast][pvalue <= max_p]) +
#       ggplot2::aes(x = .data[[x]], y = .data[["pvalue"]], colour = abs(.data[[x]])) +
#       ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "black") +
#       ggplot2::geom_point(size = 0.50) +
#       ggrepel::geom_label_repel(
#         data = ~ .x[!is.na(external_gene_name) & fdr < 0.05][order(pvalue)][1:10],
#         mapping = aes(label = .data[["external_gene_name"]]),
#         min.segment.length = unit(0, "lines"),
#         size = 1.5,
#         show.legend = FALSE,
#         na.rm = TRUE,
#         nudge_y = 0.5,
#         fontface = "italic",
#         segment.colour = "black",
#         colour = "black",
#         max.overlaps = 15
#       ) +
#       ggplot2::scale_colour_viridis_c(
#         trans = "sqrt",
#         limits = c(0, data[pvalue <= 0.05, max(abs(.SD)), .SDcols = x])
#       ) +
#       ggplot2::labs(x = label_x, y = label_y, colour = label_colour) +
#       ggplot2::theme(legend.position = "none") +
#       ggplot2::annotate(
#         geom = "rect",
#         xmin = -Inf, xmax = Inf, ymin = 1, ymax = alpha,
#         fill = "#b22222", alpha = 0.2, colour = NA
#       ) +
#       ggplot2::geom_hline(yintercept = alpha, linetype = 2, colour = "#b22222") +
#       ggplot2::scale_y_continuous(
#         trans = pval_trans(alpha = alpha, md = TRUE, colour = "#b22222"),
#         expand = ggplot2::expansion(mult = c(0, 0.2)),
#         limits = c(max_p, NA)
#       ) +
#       ggplot2::theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
#       ggplot2::theme(
#         plot.title.position = "plot",
#         plot.caption.position = "plot",
#         plot.title = ggtext::element_markdown(),
#         plot.subtitle = ggtext::element_markdown(face = "italic", size = rel(0.80)),
#         plot.caption = ggtext::element_markdown(face = "italic", size = rel(0.75)),
#         axis.title.x = ggtext::element_markdown(),
#         axis.text.x = ggtext::element_markdown(),
#         axis.title.y = ggtext::element_markdown(),
#         axis.text.y = ggtext::element_markdown(),
#         panel.grid.minor = ggplot2::element_blank(),
#         legend.title = ggtext::element_markdown()
#       )
#   })
#   names(out) <- unique(data[["contrast"]])
#   out
# }
