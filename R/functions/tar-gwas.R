#' qc_sample_sheet_gwas
#' @import data.table
qc_sample_sheet_gwas <- function(phenotype, exclusion, relatedness, ethnicity) {
  relatedness <- data.table::fread(file = relatedness)

  related_pairs <- relatedness[
    i = sub("_.*", "", IID1) != sub("_.*", "", IID2)
  ]

  best_related_samples <- sapply(
    X = unique(lapply(
      X = related_pairs[["IID1"]],
      FUN = function(iid) {
        related_pairs[IID1 %in% iid | IID2 %in% iid, sort(unique(c(IID1, IID2, iid)))]
      }
    )),
    FUN = function(x) {
      data.table::melt.data.table(
        data = related_pairs[IID1 %in% x | IID2 %in% x, .SD, .SDcols = paste0(rep(c("IID", "F_MISS"), each = 2), 1:2)],
        measure.vars = patterns("^F_MISS", "^IID"),
        value.name = c("F_MISS", "IID")
      )[
        which.min(F_MISS),
        unique(IID)
      ]
    }
  )

  bad_duplicated_samples <- relatedness[
    i = sub("_.*", "", IID1) == sub("_.*", "", IID2)
  ][
    j = data.table::melt.data.table(.SD, measure.vars = paste0("F_MISS", 1:2)),
    .SDcols = c(paste0("IID", 1:2), paste0("F_MISS", 1:2))
  ][
    j = list(
      IID = sub("_.*", "", IID1),
      best_iid = data.table::fifelse(
        test = variable[which.min(value)] == "F_MISS1",
        yes = IID1,
        no = IID2
      )
    ),
    by = c(paste0("IID", 1:2))
  ]

  exclusion <- data.table::fread(file = exclusion)[
    i = Sex_Discrepancy == 1 | Sample_Call_Rate == 1 | Sex_Missing == 1 | Heterozygosity_Check == 1,
    j = Status := "Exclude"
  ][
    i = IID %in%
      bad_duplicated_samples[j = unlist(.SD, use.names = FALSE), .SDcols = patterns("IID[0-9]+")] &
        Status == "Check",
    j = Status := NA_character_
  ]

  dt <- phenotype[j = `:=`(
    "Sample_Name" = IID,
    "vcf_id" = IID
  )][# Add related variable
    j = is_related := vcf_id %in% relatedness[
      i = sub("_.*", "", IID1) != sub("_.*", "", IID2),
      j = unique(unlist(.SD)),
      .SDcols = paste0("IID", 1:2)
    ]
  ][
    j = is_best_related := vcf_id %in% best_related_samples
  ][# keep best duplicates from genotyping
    i = bad_duplicated_samples,
    j = vcf_id := best_iid,
    on = "IID"
  ]

  merge(
    x = merge(x = dt, y = exclusion[j = list(vcf_id = IID, Status)], by = "vcf_id", all.x = TRUE),
    y = data.table::fread(file = ethnicity)[j = -c("cohort")], # Add genetics PCs
    by.x = "vcf_id",
    by.y = "iid",
    all.x = TRUE
  )[
    i = (is_related),
    j = Status := "Related"
  ][
    i = (is_related) & (is_best_related),
    j = Status := "Related (Best)"
  ][
    i = is.na(PC01),
    j = `:=`(
      vcf_id = NA_character_,
      Status = "Exclude"
    )
  ][
    i = is.na(vcf_id),
    j = is_related := NA
  ]
}

#' do_gwas
#' @import data.table
#' @import stats
#' @import utils
#' @import future.apply
#' @import rlang
do_gwas <- function(
  data,
  model,
  vcfs,
  path,
  vep = NULL,
  bin_path = list(
    bcftools = "/usr/bin/bcftools",
    plink2 = "plink2"
  )
) {
  INFO <- TEST <- P <- NULL # "no visible binding for global variable" from `data.table` syntax
  path <- normalizePath(path)
  dir.create(path = path, recursive = TRUE, mode = "0775", showWarnings = FALSE)

  plink_version <- try(
    expr = system(sprintf("%s --version", bin_path[["plink2"]]), intern = TRUE),
    silent = TRUE
  )
  if (inherits(plink_version, "try-error")) stop("Please check PLINK binary path!")

  bcftools_version <- try(
    expr = system(sprintf("%s --version", bin_path[["bcftools"]]), intern = TRUE)[1],
    silent = TRUE
  )
  if (inherits(bcftools_version, "try-error")) stop("Please check BCFTools binary path!")

  dir.create(
    path = path,
    recursive = TRUE,
    mode = "0775",
    showWarnings = FALSE
  )
  tmpdir <- file.path(tempdir(), "gwas_plink2", model[["raw_trait"]])
  dir.create(path = tmpdir, recursive = TRUE, mode = "0775")
  on.exit(unlink(tmpdir, recursive = TRUE))

  vcfs <- vcfs[sub("\\.vcf.gz", "", basename(vcfs)) %in% 1:22]

  covariates <- all.vars(as.formula(sprintf(" ~ %s", model[["covariates"]])))
  dt <- data.table::setnames(
    x = data.table::as.data.table(data)[
      j = unique(na.exclude(.SD)),
      .SDcols = c("vcf_id", model[["raw_trait"]], covariates)
    ],
    old = "vcf_id",
    new = "#IID"
  )

  if (length(sex_covariate <- grep("^sex", covariates, value = TRUE)) > 1) {
    stop(sprintf(
      'Only one column containing "sex" can exist in the model, not %s: "%s"',
      length(sex_covariate),
      paste(sex_covariate, collapse = '", "')
    ))
  }

  binary_wrongly_coded <- dt[
    j = names(which(sapply(
      X = .SD,
      FUN = function(.col) {
        data.table::uniqueN(.col) == 2 && 0 %in% unique(.col)
      }
    ))),
    .SDcols = c(model[["raw_trait"]])
  ]

  if (length(binary_wrongly_coded) > 0) {
    stop(
      "Binary traits must be coded as 1 or 2!\n",
      sprintf("Please check: \"%s\"", paste(binary_wrongly_coded, collapse = '", "'))
    )
  }

  basename_file <- file.path(tmpdir, rlang::hash(model))

  data.table::fwrite(
    x = dt[j = unique(.SD), .SDcols = "#IID"],
    file = sprintf("%s.samples", basename_file),
    sep = " ",
    col.names = FALSE
  )

  data.table::fwrite(
    x = dt[j = unique(.SD), .SDcols = c("#IID", model[["raw_trait"]])],
    file = sprintf("%s.pheno", basename_file),
    sep = " "
  )

  if (length(covariates_not_sex <- setdiff(covariates, sex_covariate)) > 0) {
    data.table::fwrite(
      x = dt[j = unique(.SD), .SDcols = c("#IID", covariates_not_sex)],
      file = sprintf("%s.cov", basename_file),
      sep = " "
    )
  }

  if (length(sex_covariate) > 0) {
    if (length(sex_levels <- unique(dt[[sex_covariate]])) == 2 & 0 %in% sex_levels) {
      stop("Sex must be coded: '1'/'M'/'m' = male, '2'/'F'/'f' = female, 'NA'/'0' = missing!")
    }
    data.table::fwrite(
      x = data.table::setnames(
        x = data.table::copy(dt)[j = unique(.SD), .SDcols = c("#IID", sex_covariate)],
        old = sex_covariate,
        new = "SEX"
      ),
      file = sprintf("%s.sex", basename_file),
      sep = " "
    )
  }

  message("Formatting VCFs and performing PLINK2 regression ...")
  if (nzchar(system.file(package = "future.apply"))) {
    gwas_lapply <- function(X, basename_file, vep_file, bin_path, FUN) {
      future.apply::future_lapply(
        X = X,
        basename_file = basename_file,
        vep_file = vep_file,
        bin_path = bin_path,
        future.globals = FALSE,
        future.packages = "data.table",
        FUN = FUN
      )
    }
  } else {
    gwas_lapply <- function(X, basename_file, vep_file, bin_path, FUN) {
      lapply(
        X = X,
        basename_file = basename_file,
        vep_file = vep_file,
        bin_path = bin_path,
        FUN = FUN
      )
    }
  }
  list_results <- gwas_lapply(
    X = vcfs,
    basename_file = basename_file,
    vep_file = vep,
    bin_path = bin_path,
    FUN = function(vcf, basename_file, vep_file, bin_path) {
      vcf_file <- sprintf("%s__%s", basename_file, basename(vcf))
      results_file <- sub("\\.vcf.gz", "", vcf_file)

      cmd <- paste(
        bin_path[["bcftools"]],
          "+fill-tags", vcf,
       "|",
        bin_path[["bcftools"]],
          "view",
          "--min-af 0.05",
          "--exclude 'INFO/INFO < 0.8'",
          "--min-alleles 2 --max-alleles 2 --types snps",
          "--force-samples",
          "--samples-file", sprintf("%s.samples", basename_file)
      )

      if (!is.null(vep_file) && file.exists(vep_file)) {
        cmd <- paste(
          cmd,
          "|",
          bin_path[["bcftools"]],
            "annotate",
            "--annotations", vep_file,
            "--header-lines", sub("_formatted.tsv.gz", ".header", vep_file),
            "--columns CHROM,POS,Gene,Symbol,rsid",
          "|",
          bin_path[["bcftools"]],
            "annotate",
            "--set-id '%INFO/rsid'",
          "|",
          bin_path[["bcftools"]],
            "annotate",
            "--set-id +'%CHROM:%POS:%REF:%ALT'",
            "--output-type z --output", vcf_file
        )
      } else {
        cmd <- paste(
          cmd,
          "|",
          bin_path[["bcftools"]],
            "annotate",
            "--set-id +'%CHROM:%POS:%REF:%ALT'",
            "--output-type z --output", vcf_file
        )
      }

      system(cmd)

      system(paste(c(
        bin_path[["plink2"]],
        "--vcf", vcf_file, "dosage=DS",
        "--mach-r2-filter",
        "--threads", threads,
        "--glm sex",
        if (file.exists(sprintf("%s.cov", basename_file))) c("--covar", sprintf("%s.cov", basename_file)) else "allow-no-covars",
        if (file.exists(sprintf("%s.samples", basename_file))) c("--keep", sprintf("%s.samples", basename_file)),
        if (file.exists(sprintf("%s.sex", basename_file))) c("--update-sex", sprintf("%s.sex", basename_file)),
        if (file.exists(sprintf("%s.pheno", basename_file))) c("--pheno", sprintf("%s.pheno", basename_file)),
        "--covar-variance-standardize",
        "--silent",
        "--out", results_file
      ), collapse = " "))

      annot <- data.table::setnames(
        x = data.table::fread(
          cmd = paste(bin_path[["bcftools"]], "view --drop-genotypes", vcf_file),
          skip = "#CHROM"
        ),
        old = function(x) sub("^#", "", x)
      )

      if (any(grepl("^INFO$", names(annot)))) {
        annot <- annot[
          j = list(
            .SD,
            data.table::rbindlist(
              l = lapply(
                X = strsplit(INFO, ";"),
                FUN = function(x) {
                  all_fields <- strsplit(x, "=")
                  out <- data.table::transpose(all_fields[sapply(all_fields, length) > 1])
                  data.table::setnames(x = data.table::setDT(do.call("rbind.data.frame", out[-1])), old = out[[1]])
                }
              ),
              use.names = TRUE,
              fill = TRUE
            )[
              j = lapply(.SD, function(x) {
                xout <- as.character(x)
                data.table::fifelse(
                  test = xout %in% c(".", "-"),
                  yes = NA_character_,
                  no = xout
                )
              })
            ]
          ),
          .SDcols = !intersect(c("INFO", "QUAL", "FILTER"), names(annot))
        ]
      }

      if (length(qual_filter_cols <- intersect(c("QUAL", "FILTER"), names(annot))) > 0) {
        annot <- annot[j = .SD, .SDcols = !c(qual_filter_cols)]
      }

      data.table::setnames(
        x = annot,
        old = function(x) sub("^\\.SD\\.\\.*", "", x)
      )

      results <- data.table::setnames(
        x = data.table::rbindlist(
          l = lapply(
            X = (function(x) `names<-`(x, sub("[^.]*\\.", "", x)))(
              list.files(
                path = dirname(results_file),
                pattern = sprintf("%s\\..*\\.glm\\..*", basename(results_file)),
                full.names = TRUE
              )
            ),
            FUN = data.table::fread
          ),
          idcol = "trait_model"
        ),
        old = function(x) sub("^#", "", x)
      )[TEST %in% "ADD" & !is.na(P), -c("TEST")]

      data.table::fwrite(
        x = data.table::merge.data.table(
          x = results,
          y = annot,
          by = c("CHROM", "POS", "ID", "REF", "ALT"), # intersect(names(results), names(annot))
          all.x = TRUE
        ),
        file = sprintf("%s.results.gz", results_file)
      )
      message(sprintf("Results written in '%.results.gz'", results_file))

      sprintf("%s.results.gz", results_file)
    }
  )

  message("Aggregating PLINK2 results ...")

  results_file <- file.path(path, "gwas.csv.gz")

  data.table::fwrite(
    x = data.table::setcolorder(
      x = data.table::rbindlist(
        l = lapply(list_results, data.table::fread),
        use.names = TRUE
      )[
        j = `:=`(
          FDR = stats::p.adjust(P, method = "BH"),
          Bonferroni = stats::p.adjust(P, method = "bonferroni"),
          covariates = covariates
        ),
        by = "trait_model"
      ][order(P)],
      neworder = c("trait_model", "covariates")
    ),
    file = results_file
  )

  message(sprintf('Writing results to "%s"!', results_file))

  results_file
}

#' plot_manhattan_gwas
#' @import ggplot2
#' @import ggtext
#' @import data.table
#' @import ggrepel
#' @import stats
plot_manhattan_gwas <- function(file, model) {
  raw_trait <- all.vars(stats::as.formula(paste0("~", model[["raw_trait"]])))

  dt <- data.table::fread(file)[
    j = gene_label_min := data.table::fifelse(
      test = P == min(P, na.rm = TRUE) &
        !is.na(Symbol) &
        !Symbol %in% c("", "NA"),
      yes = gsub(",", ";", Symbol),
      no = NA_character_
    ),
    by = "Symbol"
  ][
    i = P > 0.05,
    j = P := NA_real_
  ][
    j = file := basename(file)
  ][
    j = c("CHROM", "POS", "P", "gene_label_min")
  ][order(P)]

  if (is.numeric(dt[["CHROM"]])) {
    dt[j = "CHROM" := lapply(.SD, factor), .SDcols = "CHROM"]
  }

  if (dt[!is.na(gene_label_min), .N] > 10) {
    dt[which(!is.na(gene_label_min))[-c(1:10)], gene_label_min := NA_character_]
  }

  alpha <- 0.05 / nrow(dt)

  ggplot2::ggplot(data = dt) +
    ggplot2::aes(x = .data[["POS"]], y = .data[["P"]], colour = .data[["CHROM"]]) +
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

#' plot_pp_gwas
#' @import data.table
#' @import ggplot2
#' @import ggtext
#' @import stats
plot_pp_gwas <- function(file, model) {
  raw_trait <- all.vars(stats::as.formula(paste0("~", model[["raw_trait"]])))

  dt <- data.table::fread(file)[
    order(P)
  ][
    j = `:=`(
      "exppval" = (seq_len(.N) - 0.5) / .N,
      "labels" = paste0(
        "&lambda;<sub>gc</sub> = ",
        format(stats::median(stats::qnorm(P / 2)^2, na.rm = TRUE) / stats::qchisq(0.5, df = 1), digits = 3, nsmall = 3)
      )
    )
  ]

  alpha <- 0.05 / nrow(dt)

  ggplot2::ggplot(data = dt) +
    ggplot2::aes(
      x = .data[["exppval"]],
      y = .data[["P"]],
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

#' plot_volcano_gwas
#' @import ggplot2
#' @import ggtext
#' @import data.table
#' @import ggrepel
#' @import stats
plot_volcano_gwas <- function(file, model) {
  raw_trait <- all.vars(stats::as.formula(paste0("~", model[["raw_trait"]])))

  dt <- data.table::fread(file)[
    j = gene_label_min := data.table::fifelse(
      test = P == min(P, na.rm = TRUE) &
        !is.na(Symbol) &
        !Symbol %in% c("", "NA"),
      yes = gsub(",", ";", Symbol),
      no = NA_character_
    ),
    by = "Symbol"
  ][
    i = P > 0.05,
    j = P := NA_real_
  ][
    j = file := basename(file)
  ][
    j = .SD,
    .SDcols = patterns("^(OR|BETA|CHROM|POS|P|gene_label_min)$")
  ][order(P)]

  if (is.numeric(dt[["CHROM"]])) {
    dt[j = "CHROM" := lapply(.SD, as.character), .SDcols = "CHROM"]
  }

  if (dt[!is.na(gene_label_min), .N] > 10) {
    dt[which(!is.na(gene_label_min))[-c(1:10)], gene_label_min := NA_character_]
  }

  alpha <- 0.05 / nrow(dt)

  ggplot2::ggplot(dt) +
    ggplot2::aes(
      x = .data[[grep("^(OR|BETA)$", names(dt), value = TRUE)]],
      y = .data[["P"]],
      colour = abs(.data[[grep("^(OR|BETA)$", names(dt), value = TRUE)]])
    ) +
    ggplot2::geom_vline(xintercept = as.numeric(grepl("^OR$", names(dt))), linetype = 2) +
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
      x = c("BETA" = "Estimate", "OR" = "Odds Ratio")[grep("^(OR|BETA)$", names(dt), value = TRUE)],
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
