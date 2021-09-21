#' qc_sample_sheet_ewas
#' @import data.table
qc_sample_sheet_ewas <- function(phenotype, methylation) {
  methy_sample_sheet <- data.table::fread(file = methylation)

  merge(
    x = phenotype[j = `:=`("Sample_Name" = IID)],
    y =   methy_sample_sheet[
      j = .SD,
      .SDcols = grep("^CellT_|^Sample_|^Sentrix_|qc_sex_discrepancy|call_rate", names(methy_sample_sheet))
    ], # Add cell components and QC metrics
    by = "Sample_Name",
    all.x = TRUE
  )[
    i = (qc_sex_discrepancy | is.na(call_rate)),
    j = Status := "Exclude"
  ][
    j = if (.N > 1) .SD[!Status %in% "Exclude" & call_rate == max(call_rate, na.rm = TRUE)] else .SD,
    by = "Sample_Name"
  ]
}

#' do_ewas
#' @import data.table
#' @import stats
#' @import utils
#' @import limma
#' @import IlluminaHumanMethylationEPICanno.ilm10b5.hg38
#' @importFrom future.apply future_lapply future_apply
do_ewas <- function(
  data,
  model,
  beta_file,
  path,
  epic_annot_pkg = "IlluminaHumanMethylationEPICanno.ilm10b5.hg38"
) {
  tmpdir <- file.path(tempdir(), "ewas_limma", model[["raw_trait"]])
  dir.create(path = tmpdir, recursive = TRUE, mode = "0777")
  on.exit(unlink(tmpdir, recursive = TRUE))

  beta_matrix <- data.table::fread(file = beta_file, header = TRUE)
  data <- data[Sample_ID %in% names(beta_matrix)]
  beta_matrix <- data.table::setnames(beta_matrix, data[["Sample_ID"]], data[["Sample_Name"]])[
    j = (function(x) log2(x) - log2(1 - x))(as.matrix(.SD, "cpg_id")),
    .SDcols = c("cpg_id", data[["Sample_Name"]])
  ]

  covariates <- all.vars(stats::as.formula(sprintf(" ~ %s", model[["covariates"]])))
  if (any(grepl("^cell$", covariates))) {
    covariates <- c(
      covariates[!grepl("^cell$", covariates)],
      names(sort(data[j = colMeans(.SD, na.rm = TRUE), .SDcols = grep("^CellT_", names(data))])[-1])
    )
  }

  if (length(sex_covariate <- grep("^sex", covariates, value = TRUE)) > 1) {
    stop(sprintf(
      'Only one column containing "sex" can exist in the model, not %s: "%s"',
      length(sex_covariate),
      paste(sex_covariate, collapse = '", "')
    ))
  }

  basename_file <- sprintf("%s/limma_%s", tmpdir, model[["raw_trait"]])

  message("Performing limma regression ...")

  form <- stats::as.formula(paste0("~ ", paste(c(model[["raw_trait"]], covariates), collapse = " + ")))
  pheno_dt <- unique(data.table::setnames(
    x = data.table::as.data.table(data)[
      j = na.exclude(.SD),
      .SDcols = c("Sample_Name", all.vars(form))
    ],
    old = "Sample_Name",
    new = "#IID"
  ))

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

  limma_fit1 <- limma::lmFit(
    object = beta_matrix[, as.character(pheno_dt[["#IID"]])],
    design = model.matrix(object = form, data = pheno_dt)
  )

  if (is.factor(trait_values) & nlevels(trait_values) > 2) {
    limma_fit2 <- limma::eBayes(limma_fit1)
    limma_top1 <- future.apply::future_lapply(
      X = paste0(model[["raw_trait"]], levels(trait_values)[-1]),
      .fit = limma_fit2,
      .trait = model[["raw_trait"]],
      .ref = levels(trait_values)[1],
      .number = nrow(beta_matrix),
      future.globals = FALSE,
      future.packages = c("data.table", "limma", "stats"),
      FUN = function(.coef, .fit, .trait, .ref, .number) {
        data.table::as.data.table(limma::topTable(
          fit = .fit,
          coef = .coef,
          number = .number,
          adjust.method = "BH",
          p.value = 1,
          sort.by = "none"
        ), keep.rownames = "CpG")[
          j = `:=`(
            "se" = sqrt(.fit[["s2.post"]]) * .fit[["stdev.unscaled"]][, .coef],
            "adj.P.Val" = NULL,
            "B" = NULL,
            "fdr" = stats::p.adjust(P.Value, method = "BH"),
            "contrast" = sprintf("%s: %s Vs. %s (ref)", .trait, sub(.trait, "", .coef), .ref)
          )
        ]
      }
    )
    limma_top2 <- future.apply::future_apply(
      X = combn(levels(trait_values)[-1], 2),
      MARGIN = 2,
      .fit = limma_fit1,
      .trait = model[["raw_trait"]],
      .number = nrow(beta_matrix),
      future.globals = FALSE,
      future.packages = c("data.table", "limma", "stats"),
      FUN = function(icol, .fit, .trait, .number) {
        out <- paste0(.trait, icol)
        coef <- paste0(out[2], "-", out[1])
        contrasts_fit <- limma::makeContrasts(contrasts = coef, levels = .fit$design)
        attr(contrasts_fit, "dimnames") <- lapply(
          X = attr(contrasts_fit, "dimnames"),
          FUN = gsub, pattern = "^Intercept$", replacement = "(Intercept)"
        )
        limma_fit1b <- limma::contrasts.fit(fit = .fit, contrasts = contrasts_fit)
        limma_fit2b <- limma::eBayes(limma_fit1b)
        data.table::as.data.table(limma::topTable(
          fit = limma_fit2b,
          coef = coef,
          number = .number,
          adjust.method = "BH",
          p.value = 1,
          sort.by = "none"
        ), keep.rownames = "CpG")[
          j = `:=`(
            "se" = sqrt(limma_fit2b[["s2.post"]]) * limma_fit2b[["stdev.unscaled"]][, coef],
            "adj.P.Val" = NULL,
            "B" = NULL,
            "fdr" = stats::p.adjust(P.Value, method = "BH"),
            "contrast" = sprintf("%s: %s Vs. %s (ref)", .trait, icol[2], icol[1])
          )
        ]
      }
    )
    limma_top <- data.table::rbindlist(c(limma_top1, limma_top2))
  } else {
    limma_fit2 <- limma::eBayes(limma_fit1)
    if (is.factor(trait_values)) {
      .levels <- levels(trait_values)
      .coef <- paste0(model[["raw_trait"]], .levels[-1])
      limma_top <- data.table::as.data.table(limma::topTable(
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
          "fdr" = stats::p.adjust(P.Value, method = "BH"),
          "contrast" = paste0(model[["raw_trait"]], ": ", .levels[2], " Vs. ", .levels[1], " (ref)")
        )
      ]
    } else {
      limma_top <- data.table::as.data.table(limma::topTable(
        fit = limma_fit2,
        coef = model[["raw_trait"]],
        number = nrow(beta_matrix),
        adjust.method = "BH",
        p.value = 1,
        sort.by = "none"
      ), keep.rownames = "CpG")[
        j = `:=`(
          "se" = sqrt(limma_fit2[["s2.post"]]) * limma_fit2[["stdev.unscaled"]][, model[["raw_trait"]]],
          "adj.P.Val" = NULL,
          "B" = NULL,
          "fdr" = stats::p.adjust(P.Value, method = "BH"),
          "contrast" = model[["raw_trait"]]
        )
      ]
    }
  }

  results_file <- sprintf("%s/ewas_%s_%s.csv.gz", path, model[["raw_trait"]], model[["tar_group"]])
  dir.create(
    path = dirname(results_file),
    recursive = TRUE,
    mode = "0755",
    showWarnings = FALSE
  )

  data.table::setnames(
    x = limma_top[j = `:=`(n = nrow(pheno_dt))],
    old = c("logFC", "AveExpr", "t", "P.Value"),
    new = c("estimate", "avgmvalue_meth", "t_statistic", "pvalue"),
    skip_absent = TRUE
  )

  limma_annot <- Reduce(
    f = function(x, y) data.table::merge.data.table(x, y, by = "CpG", all.x = TRUE),
    x = suppressWarnings(list(
      limma_top,
      data.table::as.data.table(
        x =   (function(x) `names<-`(x, paste0("cpg_", names(x))))(
          x = get(utils::data("Locations", package = epic_annot_pkg))
        ),
        keep.rownames = "CpG"
      ),
      data.table::as.data.table(
        x = get(utils::data("Islands.UCSC", package = epic_annot_pkg)),
        keep.rownames = "CpG"
      ),
      data.table::as.data.table(
        x = get(utils::data("Other", package = epic_annot_pkg)),
        keep.rownames = "CpG"
      )[
        j = list(CpG, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group)
      ]
    ))
  )

  fwrite(x = limma_annot, file = results_file)

  message(sprintf('Writing results to "%s"!', results_file))

  results_file
}

#' plot_manhattan_ewas
#' @import ggplot2
#' @import ggtext
#' @import data.table
#' @import ggrepel
#' @import stats
plot_manhattan_ewas <- function(file, model) {
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
    j = c("cpg_chr", "cpg_pos", "pvalue", "gene_label_min")
  ][order(pvalue)]

  if (is.numeric(dt[["cpg_chr"]])) {
    dt[j = "cpg_chr" := lapply(.SD, as.character), .SDcols = "cpg_chr"]
  }

  if (dt[!is.na(gene_label_min), .N] > 10) {
    dt[which(!is.na(gene_label_min))[-c(1:10)], gene_label_min := NA_character_]
  }

  alpha <- 0.05 / nrow(dt)

  ggplot2::ggplot(data = dt) +
    ggplot2::aes(x = .data[["cpg_pos"]], y = .data[["pvalue"]], colour = .data[["cpg_chr"]]) +
    ggplot2::geom_point(stat = "manhattan", size = 0.60, na.rm = TRUE) +
    ggplot2::annotate(
      geom = "rect",
      xmin = -Inf, xmax = Inf, ymin = 1, ymax = alpha,
      fill = "#b22222", alpha = 0.2, colour = NA,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = alpha, linetype = 2, colour = "#b22222") +
    ggplot2::scale_x_continuous(
      breaks = 1:24,
      labels = c(1:22, "X", "Y"),
      expand = ggplot2::expansion(add = 0.25)
    ) +
    ggplot2::scale_y_continuous(
      trans = pval_trans(md = TRUE),
      expand = ggplot2::expansion(mult = c(0, 0.2))#,
      # limits = c(0.05, NA)
    ) +
    ggplot2::coord_cartesian(ylim = c(0.05, NA)) +
    ggplot2::scale_colour_manual(values = rep(scales::viridis_pal(begin = 1/4, end = 3/4)(2), 12)) +
    ggrepel::geom_label_repel(
      mapping = ggplot2::aes(label = .data[["gene_label_min"]]),
      stat = "manhattan",
      show.legend = FALSE,
      min.segment.length = 0,
      # direction = "x",
      size = 2.5,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      x = "Chromosome",
      y = "P-value",
      colour = "Chromosome",
      title = model[["pretty_trait"]],
      subtitle = toupper(paste(raw_trait, "= ", model[["covariates"]]))
    ) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggtext::element_markdown(),
      plot.subtitle = ggtext::element_markdown(face = "italic"),
      axis.text.y = ggtext::element_markdown(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.position = "none"
    )
}

#' plot_pp_ewas
#' @import data.table
#' @import ggplot2
#' @import ggtext
#' @import stats
plot_pp_ewas <- function(file, model) {
  raw_trait <- all.vars(stats::as.formula(paste0("~", model[["raw_trait"]])))

  dt <- data.table::fread(file)[
    order(pvalue)
  ][
    j = `:=`(
      "exppval" = (1:.N - 0.5) / .N,
      "labels" = paste0(
        "&lambda;<sub>gc</sub> = ",
        format(stats::median(stats::qnorm(pvalue / 2)^2, na.rm = TRUE) / stats::qchisq(0.5, df = 1), digits = 3, nsmall = 3)
      )
    )
  ]

  alpha <- 0.05 / nrow(dt)

  ggplot2::ggplot(data = dt) +
    ggplot2::aes(
      x = .data[["exppval"]],
      y = .data[["pvalue"]],
      colour = .data[["labels"]],
      shape = .data[["labels"]]
    ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "black", linetype = 2) +
    ggplot2::geom_point(size = 0.60) +
    ggplot2::geom_hline(yintercept = alpha, linetype = 2, colour = "#b22222") +
    ggplot2::scale_x_continuous(
      trans = pval_trans(md = TRUE),
      expand = ggplot2::expansion(c(0, 0.2)),
      limits = c(1, NA)
    ) +
    ggplot2::scale_y_continuous(
      trans = pval_trans(md = TRUE),
      expand = ggplot2::expansion(c(0, 0.2)),
      limits = c(1, NA)
    ) +
    ggplot2::scale_colour_viridis_d(begin = 0.5, end = 0.5) +
    ggplot2::scale_shape_discrete(solid = TRUE) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 2))) +
    ggplot2::labs(
      x = "Expected P-value",
      y = "Observed P-value",
      colour = NULL,
      shape = NULL,
      title = model[["pretty_trait"]],
      subtitle = toupper(paste(raw_trait, "= ", model[["covariates"]]))
    ) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggtext::element_markdown(),
      plot.subtitle = ggtext::element_markdown(face = "italic"),
      axis.text.x = ggtext::element_markdown(),
      axis.text.y = ggtext::element_markdown(),
      legend.position = c(0.99, 0.01),
      legend.justification = c("right", "bottom"),
      legend.box.just = "right",
      legend.text = ggtext::element_markdown(),
      legend.margin = ggplot2::margin(1.5, 1.5, 1.5, 1.5),
      legend.spacing.x = ggplot2::unit(1.5, "pt"),
      legend.spacing.y = ggplot2::unit(1.5, "pt")
    )
}

#' plot_volcano_ewas
#' @import ggplot2
#' @import ggtext
#' @import data.table
#' @import ggrepel
#' @import stats
plot_volcano_ewas <- function(file, model) {
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
