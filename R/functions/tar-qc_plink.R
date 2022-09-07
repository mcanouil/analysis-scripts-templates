#' create_ga_directory
create_ga_directory <- function(path) {
  unlink(x = file.path(path, c("tmp", "plink_qc", "vcf_qc")), force = TRUE, recursive = TRUE)
  invisible(
    sapply(
      X = file.path(path, c("tmp", "plink_qc", "vcf_qc")),
      FUN = dir.create,
      showWarnings = FALSE,
      recursive = TRUE,
      mode = "0777"
    )
  )
  file.path(path, "tmp")
}

#' download_plink
#' @import utils
download_plink <- function(
  url = "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip",
  output_directory = tempdir()
) {
  out_file <- list.files(path = output_directory, pattern = "^plink$", full.names = TRUE)

  if (length(out_file) == 0 || !file.exists(out_file)) {
    utils::download.file(url = url, destfile = file.path(output_directory, basename(url)))
    files_extracted <- utils::unzip(
      zipfile = file.path(output_directory, basename(url)),
      exdir = output_directory
    )
    on.exit(unlink(
      x = c(
        file.path(output_directory, basename(url)),
        setdiff(files_extracted, grep("\\plink$", files_extracted, value = TRUE))
      ),
      force = TRUE
    ))
  }

  out_file <- list.files(path = output_directory, pattern = "^plink$", full.names = TRUE)

  Sys.chmod(out_file, mode = "0775")

  out_file
}

#' download_perl_preimputation_check
#' @import utils
download_perl_preimputation_check <- function(
  url = "https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.13-NoReadKey.zip",
  output_directory = tempdir()
) {
  out_file <- sub("\\.zip$", ".pl", file.path(output_directory, basename(url)))

  if (!file.exists(out_file)) {
    utils::download.file(url = url, destfile = file.path(output_directory, basename(url)))
    files_extracted <- utils::unzip(
      zipfile = file.path(output_directory, basename(url)),
      exdir = output_directory
    )
    perl_script <- normalizePath(grep("\\.pl$", files_extracted, value = TRUE))
    file.rename(perl_script, sub("\\.zip$", ".pl", file.path(output_directory, basename(url))))
  }

  if (file.exists(file.path(output_directory, basename(url)))) {
    unlink(x = file.path(output_directory, basename(url)), force = TRUE)
  }
  if (file.exists(file.path(output_directory, "LICENSE.txt"))) {
    unlink(x = file.path(output_directory, "LICENSE.txt"), force = TRUE)
  }

  Sys.chmod(out_file, mode = "0775")

  out_file
}

#' make_ga_bed
#' @import data.table
make_ga_bed <- function(input, output, plink) {
  if (length(list.files(dirname(input), pattern = ".bed")) != 0) {
    cmd <- paste(plink,
      "--bfile", input,
      "--snps-only",
      "--make-bed",
      "--out", file.path(output, "data")
    )
  } else {
    cmd <- paste(plink,
      "--file", input,
      "--snps-only",
      "--make-bed",
      "--out", file.path(output, "data")
    )
  }

  system(command = cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  data.table::fwrite(
    x = data.table::fread(file.path(output, "data.fam"))[,
      V2 := {
        x <- as.character(seq_len(.N))
        data.table::fifelse(x == "1", as.character(V2), paste(V2, x, sep = "_"))
      },
      by = c("V1", "V2")
    ],
    file = file.path(output, "data.fam"),
    col.names = FALSE,
    sep = "\t"
  )
  out_files <- list.files(
    path = output,
    pattern = paste(paste0("data.", c("bed", "bim", "fam")), collapse = "|"),
    full.names = TRUE
  )

  unlink(
    x = setdiff(
      list.files(path = output, pattern = "data\\.", full.names = TRUE),
      out_files
    ),
    force = TRUE
  )

  out_files
}

#' count_duplicated_samples
count_duplicated_samples <- function(data) {
  sum(table(gsub("[_-][^-_]*$", "", data[["IID"]])) > 1)
}

#' read_fam
#' @import data.table
read_fam <- function(path, project) {
  data.table::fread(
    file = path,
    col.names = c("FID", "IID", "father_id", "mother_id", "sex", "phenotype"),
    colClasses = c("character", "character", "character", "character", "numeric", "numeric")
  )[
    j = `:=`(
      "cohort" = project,
      "sex_fct" = factor(sex, levels = c(1, 2, 0), labels = c("Male", "Female", "Unspecified"))
    )
  ]
}

#' plot_callrate
#' @import data.table
#' @import ggplot2
#' @import scales
#' @import ggrepel
plot_callrate <- function(data, callrate, max_labels, type) {
  ggplot2::ggplot(data = data) +
    ggplot2::aes(x = rev(seq_along(F_MISS)), y = 1 - F_MISS) +
    ggplot2::geom_point(
      colour = scales::viridis_pal(begin = 0.5, end = 0.5)(1),
      shape = 1,
      na.rm = TRUE,
      size = 3
    ) +
    ggrepel::geom_label_repel(
      data = ~ .x[j = labs := if (sum(!is.na(labs)) > max_labels) NA else labs],
      mapping = ggplot2::aes(label = labs),
      min.segment.length = ggplot2::unit(0, "lines"),
      force = 10,
      fill = "white",
      colour  = "#b22222",
      segment.colour = "#b22222",
      size = 5,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      mapping = ggplot2::aes(yintercept = callrate),
      colour = "#b22222",
      linetype = 2
    ) +
    ggplot2::scale_x_continuous(labels = scales::comma_format(accuracy = 1), trans = "log10") +
    ggplot2::scale_y_continuous(
      breaks = function(x) {
        unique(round(c(scales::breaks_extended()(c(x, callrate)), callrate), digits = 4))
      },
      labels = function(x) {
        ifelse(
          x == callrate,
          paste0(
            "<b style='color:#b22222;'>",
            scales::percent_format(accuracy = 0.01, suffix = " %")(x),
            "</b>"
          ),
          scales::percent_format(accuracy = 0.01, suffix = " %")(x)
        )
      },
      limits = c(NA, 1)
    ) +
    switch(EXPR = type,
      "snp" = {
        ggplot2::labs(
          x = "Number of Variants",
          y = "Call Rate",
          title = "Variants Call Rate",
          caption = "Variants with a call rate greater or equal to 99.9 % are not shown."
        )
      },
      "sample" = {
        ggplot2::labs(
          x = "Number of Samples",
          y = "Call Rate",
          title = "Sample Call Rate"
        )
      }
    )
}

#' compute_callrate_ind
#' @import data.table
compute_callrate_ind <- function(bfile, callrate, plink) {
  temp_file <- tempfile(pattern = "crind")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))

  system(paste(plink,
    "--bfile", bfile,
    "--missing",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  data.table::fread(
    file = sprintf("%s.imiss", temp_file),
    colClasses = c("FID" = "character", "IID" = "character")
  )[
    j = labs := data.table::fifelse(F_MISS > 1 - callrate, IID, NA_character_)
  ][order(F_MISS)]
}

#' compute_snps_maf_leq_geq
#' @import data.table
compute_snps_maf_leq_geq <- function(bfile, maf, callrate, plink) {
  temp_file <- tempfile(pattern = "snps_maf")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))
  system(paste(plink,
    "--bfile", bfile,
    "--freq",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  data.table::fread(sprintf("%s.frq", temp_file))[!is.na(MAF)][
    j = grp := factor(MAF >= maf, levels = c(FALSE, TRUE), labels = c("leq", "geq"))
  ][
    j = (function(bfile, snps, what, plink) {
      temp_file <- tempfile(pattern = sprintf("snps_maf_%s", unique(what)))
      on.exit(unlink(
        x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
        force = TRUE
      ))
      data.table::fwrite(
        x = list(snps),
        file = sprintf("%s.txt", temp_file),
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE
      )
      if (length(snps) > 0) {
        system(paste(plink,
          "--bfile", bfile,
          "--missing",
          "--extract", sprintf("%s.txt", temp_file),
          "--out", temp_file
        ), ignore.stdout = TRUE, ignore.stderr = TRUE)

        data.table::fread(
          file = sprintf("%s.imiss", temp_file),
          colClasses = c("FID" = "character", "IID" = "character")
        )[
          j = labs := data.table::fifelse(F_MISS > 1 - callrate, IID, NA_character_)
        ][order(F_MISS)]
      }
    })(bfile = bfile, snps = SNP, what = grp, plink = plink),
    by = "grp"
  ]
}

#' check_genotypic_sex
#' @import data.table
check_genotypic_sex <- function(bfile, callrate_data, fam_data, sex_threshold, plink) {
  temp_file <- tempfile(pattern = "check_sex")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))
  system(paste(plink,
    "--bfile", bfile,
    "--check-sex", min(sex_threshold), max(sex_threshold),
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  merge(
    x = merge(
      x = data.table::fread(
        file = sprintf("%s.sexcheck", temp_file),
        colClasses = c("FID" = "character", "IID" = "character")
      ),
      y = callrate_data,
      by = c("FID", "IID"),
      all = TRUE
    ),
    y = fam_data,
    by = c("FID", "IID"),
    all.x = TRUE
  )[
    j = sex_discrepancy := sex != SNPSEX
  ][
    j = STATUS := data.table::fifelse(
      test = sex_discrepancy | is.na(sex_discrepancy) | `F` >= min(sex_threshold) & `F` < max(sex_threshold),
      yes = "PROBLEM",
      no = STATUS
    )
  ][
    j = labs := data.table::fifelse(
      test = STATUS == "PROBLEM" | !is.na(labs),
      yes = IID,
      no = NA_character_
    )
  ]
}

#' plot_check_genotypic_sex
#' @import data.table
#' @import ggplot2
#' @import scales
#' @import ggrepel
plot_check_genotypic_sex <- function(data, callrate, sex_threshold, max_labels) {
  ggplot2::ggplot(data = data[!is.na(`F`), ]) +
    ggplot2::aes(x = `F`, y = 1 - F_MISS) +
    ggplot2::geom_rect(
      mapping = ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = callrate),
      fill = "grey95",
      colour = "transparent"
    ) +
    ggplot2::annotate(
      geom = "rect",
      xmin = -Inf, xmax = -0.2, ymin = -Inf, ymax = Inf,
      fill = "#b22222",
      alpha = 0.10,
      colour = "transparent"
    ) +
    ggplot2::annotate(
      geom = "rect",
      xmin = -0.2, xmax = min(sex_threshold), ymin = -Inf, ymax = Inf,
      fill = "#b22222",
      alpha = 0.2,
      colour = "transparent"
    ) +
    ggplot2::annotate(
      geom = "rect",
      xmin = max(sex_threshold), xmax = 1, ymin = -Inf, ymax = Inf,
      fill = "#2222b2",
      alpha = 0.2,
      colour = "transparent"
    ) +
    ggplot2::geom_hline(yintercept = callrate, colour = "#b22222", linetype = 2) +
    ggplot2::geom_point(
      data = ~ .x[!(sex_discrepancy)],
      na.rm = TRUE,
      shape = 1,
      colour = scales::viridis_pal(begin = 0.5, end = 0.5)(1),
      size = 3
    ) +
    ggplot2::geom_point(
      data = ~ .x[(sex_discrepancy)],
      na.rm = TRUE,
      shape = 3,
      colour  = "#b22222",
      size = 3
    ) +
    ggrepel::geom_label_repel(
      data = ~ .x[(sex_discrepancy)][, labs := if (sum(sex_discrepancy) > max_labels) NA else labs],
      mapping = ggplot2::aes(label = labs),
      min.segment.length = ggplot2::unit(0, "lines"),
      force = 10,
      fill = "white",
      colour  = "#b22222",
      segment.colour = "#b22222",
      size = 5,
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, sex_threshold, 1),
      labels = scales::percent_format(accuracy = 0.01, suffix = " %")
    ) +
    ggplot2::scale_y_continuous(
      breaks = function(x) {
        unique(round(c(scales::breaks_extended()(c(x, callrate)), callrate), digits = 4))
      },
      labels = function(x) {
        ifelse(
          x == callrate,
          paste0(
            "<b style='color:#b22222;'>",
            scales::percent_format(accuracy = 0.01, suffix = " %")(x),
            "</b>"
          ),
          scales::percent_format(accuracy = 0.01, suffix = " %")(x)
        )
      },
      limits = c(NA, 1)
    ) +
    ggplot2::labs(
      x = "Homozygosity Rate",
      y = "Call Rate",
      title = "Homozygosity Rate Based on X-chromosome Variants",
      subtitle = paste(
        "The expected homozygosity rate range:",
        "<b style = 'color:#b22222;'>for female</b>, and",
        "<b style = 'color:#2222b2;'>for male</b>.<br>",
        "Samples for which sex information is discordant or missing are denoted by \"<b style = 'color:#b22222;'>+</b>\"."
      )
    ) +
    ggplot2::theme(legend.position = "none")
}

#' compute_snps_het_all
#' @import data.table
compute_snps_het_all <- function(...) {
  data.table::rbindlist(list(...), use.names = TRUE)
}

#' compute_snps_het
#' @import data.table
compute_snps_het <- function(
  bfile,
  heterozygosity_threshold,
  callrate,
  callrate_data,
  plink,
  snps = NULL,
  what = NULL
) {
  if (is.null(snps) && length(snps) == 0 && is.null(what)) what <- "all"

  temp_file <- tempfile(pattern = sprintf("snps_het_%s", unique(what)))
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))

  if (!is.null(snps) && length(snps) > 0) {
    data.table::fwrite(
      x = list(snps),
      file = sprintf("%s.txt", temp_file),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
    system(paste(plink,
      "--bfile", bfile,
      "--het",
      "--extract", sprintf("%s.txt", temp_file),
      "--out", temp_file
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  } else {
    system(paste(plink,
      "--bfile", bfile,
      "--het",
      "--out", temp_file
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  }

  out <- merge(
    x = data.table::fread(
      file = sprintf("%s.het", temp_file),
      colClasses = c("FID" = "character", "IID" = "character"),
      check.names = TRUE
    ),
    y = callrate_data,
    by = c("FID", "IID"),
    all = TRUE
  )[
    j = hrate := (N.NM. - O.HOM.) / N.NM.,
    by = F_MISS < 1 - callrate
  ][
    j = colour := sqrt(
      (hrate - mean(hrate, na.rm = TRUE))^2 +
        (F_MISS - mean(F_MISS, na.rm = TRUE))^2
    ),
    by = F_MISS < 1 - callrate
  ][
    j = labs := data.table::fifelse(
      test = F_MISS < 1 - callrate & (
        hrate > mean(hrate, na.rm = TRUE) +
          heterozygosity_threshold * stats::sd(hrate, na.rm = TRUE) |
        hrate < mean(hrate, na.rm = TRUE) -
          heterozygosity_threshold * stats::sd(hrate, na.rm = TRUE)
      ),
      yes = IID,
      no = NA_character_
    ),
    by = F_MISS < 1 - callrate
  ]

  if (what == "all") out[j = grp := what]

  out
}

#' compute_snps_het_leq_geq
#' @import data.table
compute_snps_het_leq_geq <- function(
  bfile,
  heterozygosity_threshold,
  maf,
  callrate,
  callrate_data,
  plink
) {
  temp_file <- tempfile(pattern = "snps_het_lg")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))
  system(paste(plink,
    "--bfile", bfile,
    "--freq",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  data.table::fread(sprintf("%s.frq", temp_file))[!is.na(MAF)][
    j = grp := factor(MAF >= maf, levels = c(FALSE, TRUE), labels = c("leq", "geq"))
  ][
    j = compute_snps_het(
      bfile = bfile,
      heterozygosity_threshold = heterozygosity_threshold,
      callrate = callrate,
      callrate_data = callrate_data,
      plink = plink,
      snps = SNP,
      what = grp
    ),
    by = "grp"
  ]
}

#' plot_het_ind
#' @import data.table
#' @import ggplot2
#' @import scales
#' @import patchwork
plot_het_ind <- function(data, heterozygosity_threshold, callrate, maf, max_labels) {
  geoms_list <- list(
    ggplot2::aes(x = 1 - F_MISS, y = hrate),
    ggplot2::geom_rect(
      mapping = ggplot2::aes(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = callrate),
      fill = "grey95",
      colour = "transparent"
    ),
    ggplot2::geom_point(
      mapping = ggplot2::aes(colour = colour),
      shape = 1,
      na.rm = TRUE,
      size = 3
    ),
    ggplot2::geom_vline(
      mapping = ggplot2::aes(xintercept = callrate),
      colour = "#b22222",
      linetype = 2
    ),
    ggplot2::geom_hline(
      data = ~ .x[
        F_MISS < 1 - callrate,
        list(
          yintercept = mean(hrate, na.rm = TRUE) -
            heterozygosity_threshold * stats::sd(hrate, na.rm = TRUE)
        )
      ],
      mapping = ggplot2::aes(yintercept = yintercept),
      colour = "#b22222",
      linetype = 2
    ),
    ggplot2::geom_hline(
      data = ~ .x[
        F_MISS < 1 - callrate,
        list(
          yintercept = mean(hrate, na.rm = TRUE) +
            heterozygosity_threshold * stats::sd(hrate, na.rm = TRUE)
        )
      ],
      mapping = ggplot2::aes(yintercept = yintercept),
      colour = "#b22222",
      linetype = 2
    ),
    ggrepel::geom_label_repel(
      data = ~ .x[!is.na(labs)],
      mapping = ggplot2::aes(label = labs),
      min.segment.length = ggplot2::unit(0, "lines"),
      force = 10,
      fill = "white",
      colour  = "#b22222",
      segment.colour = "#b22222",
      size = 5,
      na.rm = TRUE
    ),
    ggplot2::scale_colour_viridis_c(begin = 0.2, end = 0.8, trans = "sqrt"),
    ggplot2::scale_x_continuous(
      breaks = function(x) {
        unique(round(c(scales::breaks_extended()(c(x, callrate)), callrate), digits = 4))
      },
      labels = function(x) {
        ifelse(
          x == callrate,
          paste0(
            "<b style='color:#b22222;'>",
            scales::percent_format(accuracy = 0.01, suffix = " %")(x),
            "</b>"
          ),
          scales::percent_format(accuracy = 0.01, suffix = " %")(x)
        )
      },
      # guide = ggplot2::guide_axis(angle = 45),
      limits = c(NA, 1)
    ),
    ggplot2::labs(x = "Call Rate", y = "Heterozygosity Rate"),
    ggplot2::theme(legend.position = "none")
  )

  patchwork::wrap_plots(
    ggplot2::ggplot(
      data = data[grp %in% "all"][j = labs := if (sum(!is.na(labs)) > max_labels) NA else labs]
    ) + geoms_list,
    ggplot2::ggplot(
      data = data[grp %in% "geq"][j = labs := if (sum(!is.na(labs)) > max_labels) NA else labs]
    ) + geoms_list,
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = "Heterozygosity Rate",
      tag_levels = "A",
      subtitle = paste0(
        "The horizontal red lines (\"<b style = 'color:#b22222;'>---</b>\") denotes the ", heterozygosity_threshold,
        " times the standard deviations (SD) from the mean heterozygosity rate,<br>",
        "based on: <b>A</b>) all variants, <b>B</b>) variants with a MAF &ge; ", maf, "."
      )
    )
}

#' compute_relatedness
#' @import data.table
compute_relatedness <- function(bfile, maf, ibd, plink) {
  temp_file <- tempfile(pattern = "snps_relatedness")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))

  system(paste(plink,
    "--bfile", bfile,
    "--autosome",
    "--set-hh-missing",
    "--maf", maf,
    "--r2",
    "--ld-window-r2 0.2",
    "--ld-window-kb 50",
    "--make-bed",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  system(paste(plink,
  "--bfile", temp_file,
  "--freq",
  "--out ", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  system(paste(plink,
    "--bfile", temp_file,
    "--read-freq", sprintf("%s.frq", temp_file),
    "--genome",
    "--min", ibd,
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  unique(
    data.table::fread(
      file = sprintf("%s.genome", temp_file),
      colClasses = c(
        "FID1" = "character", "IID1" = "character",
        "FID2" = "character", "IID2" = "character",
        "PI_HAT" = "numeric"
      )
    )
  )
}

#' compute_pca_ethnicity
#' @import data.table
#' @import future.apply
#' @import flashpcaR
compute_pca_ethnicity <- function(bfile, maf, hwe, ref1kg_genotypes, plink) {
  temp_file <- tempfile(pattern = "snps_ethnicity")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))

  # Target dataset
  system(paste(plink,
    "--bfile", bfile,
    "--autosome",
    "--snps-only",
    "--maf", maf,
    "--hwe", hwe,
    "--geno 0.1",
    "--make-bed",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  system(paste(plink,
    "--bfile", temp_file,
    "--freq",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  data.table::fwrite(
    x = data.table::fread(sprintf("%s.frq", temp_file))[j = list(SNP, P = 1 - MAF)],
    file = sprintf("%s.clump_maf", temp_file),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  system(paste(plink,
    "--bfile", temp_file,
    "--clump", sprintf("%s.clump_maf", temp_file),
    "--clump-p1 1",
    "--clump-p2 1",
    "--clump-r2 0.2",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  # Reference dataset
  clumped_ref_files <- future.apply::future_sapply(
    X = 1:22,
    FUN = function(ichr, ref1kg_genotypes, temp_file, plink) {
      system(paste(plink,
        "--bfile", paste0(ref1kg_genotypes, "Chr", ichr),
        "--extract", sprintf("%s.clumped", temp_file),
        "--make-bed",
        "--out", sprintf("%s_ref%s", temp_file, ichr)
      ), ignore.stdout = TRUE, ignore.stderr = TRUE)
      sprintf("%s_ref%s", temp_file, ichr)
    },
    ref1kg_genotypes = ref1kg_genotypes,
    temp_file = temp_file,
    plink = plink
  )

  data.table::fwrite(
    x = list(clumped_ref_files),
    file = sprintf("%s_ref.bmerge_list", temp_file),
    col.names = FALSE
  )

  system(paste(plink,
    "--merge-list", sprintf("%s_ref.bmerge_list", temp_file),
    "--snps-only",
    "--make-bed",
    "--out", sprintf("%s_reference", temp_file)
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  # Combine target and reference
  system(paste(plink,
    "--bfile", bfile,
    "--extract", sprintf("%s_reference.bim", temp_file),
    "--make-bed",
    "--out", sprintf("%s_target", temp_file)
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Flip to fix possible strand inconsistency
  system(paste(plink,
    "--bfile", sprintf("%s_reference", temp_file),
    "--bmerge", sprintf("%s_target", temp_file),
    "--make-bed",
    "--out", sprintf("%s_target_reference", temp_file)
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (file.exists(sprintf("%s_target_reference-merge.missnp", temp_file))) {
    system(paste(plink,
      "--bfile", sprintf("%s_target", temp_file),
      "--flip", sprintf("%s_target_reference-merge.missnp", temp_file),
      "--make-bed",
      "--out", sprintf("%s_target_flip", temp_file)
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
    unlink(sprintf("%s_target_reference-merge.missnp", temp_file))
    system(paste(plink,
      "--bfile", sprintf("%s_reference", temp_file),
      "--bmerge", sprintf("%s_target_flip", temp_file),
      "--make-bed",
      "--out", sprintf("%s_target_reference", temp_file)
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  }

  ## Exclude remaining strand inconsistency or multiallelic alleles
  if (file.exists(sprintf("%s_target_reference-merge.missnp", temp_file))) {
    if (file.exists(sprintf("%s_target_flip.bim", temp_file))) {
      input_file <- sprintf("%s_target_flip", temp_file)
    } else {
      input_file <- sprintf("%s_target", temp_file)
    }
    system(paste(plink,
      "--bfile", input_file,
      "--exclude", sprintf("%s_target_reference-merge.missnp", temp_file),
      "--make-bed",
      "--out", sprintf("%s_target_final", temp_file)
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste(plink,
      "--bfile", sprintf("%s_reference", temp_file),
      "--exclude", sprintf("%s_target_reference-merge.missnp", temp_file),
      "--make-bed",
      "--out", sprintf("%s_reference_final", temp_file)
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
    unlink(sprintf("%s_target_reference-merge.missnp", temp_file))
    system(paste(plink,
      "--bfile", sprintf("%s_reference_final", temp_file),
      "--bmerge", sprintf("%s_target_final", temp_file),
      "--make-bed",
      "--out", sprintf("%s_target_reference", temp_file)
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  }

  # Remove palindromic variant
  data.table::fwrite(
    x = data.table::fread(sprintf("%s_target_reference.bim", temp_file))[
      i = paste0(V5, V6) %in% c("GC", "CG", "AT", "TA"),
      j = list(V2)
    ],
    file = sprintf("%s_target_final.palindromic", temp_file),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  system(paste(plink,
    "--bfile", sprintf("%s_target_reference", temp_file),
    "--exclude", sprintf("%s_target_final.palindromic", temp_file),
    "--make-bed",
    "--out", sprintf("%s_target_reference_final", temp_file)
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  fid_iid <- data.table::fread(
    file = sprintf("%s_target_reference_final.fam", temp_file),
    colClasses = c("V1" = "character", "V2" = "character"),
    header = FALSE
  )

  flashpcaR::flashpca(
    X = sprintf("%s_target_reference_final", temp_file),
    ndim = 10,
    do_loadings = TRUE
  )
}

#' compute_ethnicity
#' @import data.table
tidy_pca_ethnicity <- function(data, ref1kg_panel) {
  ref_pop_table <- data.table::fread(
    file = ref1kg_panel,
    colClasses = "character",
    header = FALSE,
    select = 1:3,
    col.names = c("iid", "pop", "super_pop")
  )[j = "cohort" := "1,000 Genomes"]

  pca_dfxy <- merge(
    x = (function(vec) {
      vec_dt <- data.table::as.data.table(vec, keep.rownames = "iid")
      data.table::setnames(
        x = vec_dt,
        old = setdiff(names(vec_dt), "iid"),
        new = sprintf("PC%02d", as.numeric(gsub("V", "", setdiff(names(vec_dt), "iid"))))
      )[
        j = c("FID", "IID") := data.table::tstrsplit(x = iid, ":")
      ][
        j = "iid" := data.table::fifelse(
          test = IID == FID | !any(duplicated(IID)),
          yes = IID,
          no = paste(FID, IID, sep = "_")
        )
      ][
        j = -c("FID", "IID")
      ]
    })(data),
    y = ref_pop_table,
    by = "iid",
    all = TRUE
  )[
    i = !cohort %in% "1,000 Genomes",
    j = (c("pop", "super_pop", "cohort")) := as.list(rep("Target", 3))
  ][
    j = (c("pop", "super_pop", "cohort")) := lapply(
      X = .SD,
      FUN = function(x) factor(x, levels = c("Target", sort(setdiff(x, "Target"))))
    ),
    .SDcols = c("pop", "super_pop", "cohort")
  ]

  pop_centre <- pca_dfxy[
    i = cohort %in% "1,000 Genomes",
    j = lapply(.SD, mean),
    .SDcols = sprintf("PC%02d", 1:2),
    by = "pop"
  ]

  pca_dfxy[
    j = "pop_closest" := (function(.x, .y) {
      as.character(
        pop_centre[
          j = list("dist" = sqrt((.x - PC01)^2 + (.y - PC02)^2)),
          by = "pop"
        ][
          j = "which_closest" := dist == min(dist)
        ][(which_closest)][["pop"]]
      )
    })(PC01, PC02),
    by = "iid"
  ]

  setcolorder(
    x = merge(
      x = pca_dfxy,
      y = unique(ref_pop_table[, list("pop_closest" = pop, "super_pop_closest" = super_pop)]),
      by = "pop_closest",
      all.x = TRUE
    )[order(cohort)],
    neworder = c(
      "iid", "cohort",
      "super_pop_closest", "pop_closest",
      "super_pop", "pop"
    )
  )
}

#' compute_ethnicity
#' @import data.table
#' @import ggplot2
#' @import ggforce
#' @import concaveman
#' @import ggtext
#' @import patchwork
#' @import scales
plot_pca_ethnicty <- function(data, pve, loadings) {
  pca_contrib <- (function(pve) {
    pve_dt <- sprintf(
      fmt = "PC%02d (%s %%)",
      seq_along(pve),
      format(pve * 100, digits = 2, nsmall = 2, trim = TRUE)
    )
    names(pve_dt) <- sprintf("PC%02d", seq_along(pve))
    pve_dt
  })(pve)

  p1 <- ggplot2::ggplot(
    data = data.table::setnames(
      x = data[
        j = .SD[
          pop %in% .SD[!cohort %in% "1,000 Genomes", unique(pop_closest)] |
            !cohort %in% "1,000 Genomes"
        ]
      ],
      old = names(pca_contrib),
      new = pca_contrib
    )
  ) +
    ggplot2::aes(
      x = .data[[pca_contrib[["PC01"]]]],
      y = .data[[pca_contrib[["PC02"]]]],
      colour = pop, fill = pop, label = pop
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 0.5, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, size = 0.5, na.rm = TRUE) +
    ggforce::geom_mark_hull(
      data = ~ .x[cohort %in% "1,000 Genomes"],
      concavity = 2,
      expand = ggplot2::unit(1, "mm"),
      radius = ggplot2::unit(1, "mm"),
      con.cap = ggplot2::unit(0, "mm"),
      con.arrow = ggplot2::arrow(angle = 45, length = ggplot2::unit(2, "mm")),
      label.colour = "grey50",
      con.colour = "grey50",
      label.fontsize = 12,
      label.buffer = ggplot2::unit(2.5, "mm")
    ) +
    ggplot2::geom_point(
      data = ~ .x[!cohort %in% "1,000 Genomes"],
      fill = scales::viridis_pal(begin = 0.5, end = 0.5)(1),
      colour = "white",
      shape = 21,
      size = 3
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::scale_colour_viridis_d(na.translate = FALSE, drop = FALSE, begin = 0.10, end = 0.90) +
    ggplot2::scale_fill_viridis_d(na.translate = FALSE, drop = FALSE, begin = 0.10, end = 0.90) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position = "none"
    )

  p2 <- p1 + ggplot2::aes(colour = super_pop, fill = super_pop, label = super_pop)

  patchwork::wrap_plots(p1, p2, ncol = 2) +
    patchwork::plot_annotation(
      title = "Ethnicity Inference Based On 1,000 Genomes Project Data",
      subtitle = paste0(
        "Principal Component Analysis using ",
        format(nrow(loadings), big.mark = ",", digits = 1L, nsmall = 0L),
        " SNPs, with <b>A</b>) population level, and <b>B</b>) super population level."
      ),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.subtitle = ggtext::element_markdown(face = "italic", size = ggplot2::rel(0.80))
      )
    )
}

#' compute_samples_to_exclude
#' @import data.table
compute_samples_to_exclude <- function(
  callrate_data,
  sexcheck_data,
  heterozygosity_data,
  relatedness_data
) {
  all_exclusion <- rbind(
    callrate_data[
      i = !is.na(labs),
      j = list(FID, IID, QC = "Sample_Call_Rate")
    ],
    sexcheck_data[
      i = STATUS %in% "PROBLEM"
    ][
      i = (sex_discrepancy & !(PEDSEX %in% 0 | SNPSEX %in% 0)),
      j = list(FID, IID, QC = "Sex_Discrepancy")
    ],
    sexcheck_data[
      i = STATUS %in% "PROBLEM"
    ][
      i = is.na(sex_discrepancy) | PEDSEX %in% 0 | SNPSEX %in% 0,
      j = list(FID, IID, QC = "Sex_Missing")
    ],
    sexcheck_data[
      i = STATUS %in% "PROBLEM"
    ][
      i = !(sex_discrepancy),
      j = list(FID, IID, QC = "Sex_Frange")
    ],
    heterozygosity_data[
      i = !is.na(labs) & grp %in% "all",
      j = list(FID, IID, QC = "Heterozygosity_Check")
    ],
    heterozygosity_data[
      i = !is.na(labs) & grp %in% "geq",
      j = list(FID, IID, QC = "Heterozygosity_Check_Upper_MAF")
    ],
    unique(melt(
      data = relatedness_data[
        i = !paste(FID1, IID1, sep = "_") %in%
            callrate_data[!is.na(labs), paste(FID, IID, sep = "_")] &
          !paste(FID2, IID2, sep = "_") %in%
            callrate_data[!is.na(labs), paste(FID, IID, sep = "_")],
        j = c("FID1", "IID1", "FID2", "IID2")
      ],
      measure.vars = patterns(FID = "FID.*", IID = "IID.*")
    )[j = list(FID, IID, QC = "Relatedness")])
  )[
    j = `:=`(
      "QC" = factor(
        x = QC,
        levels = c(
          "Sample_Call_Rate",
          "Sex_Discrepancy",
          "Sex_Missing",
          "Sex_Frange",
          "Heterozygosity_Check",
          "Heterozygosity_Check_Upper_MAF",
          "Relatedness"
        )
      ),
      "value" = 1L
    )
  ]

  all_exclusion_wide <- data.table::dcast(
    data = all_exclusion,
    formula = FID + IID ~ QC,
    value.var = "value",
    fill = 0L
  )

  missing_qc_columns <- setdiff(levels(all_exclusion[["QC"]]), colnames(all_exclusion_wide))
  if (length(missing_qc_columns) > 0) {
    for (icol in missing_qc_columns) {
      all_exclusion_wide[j = tmp_col := 0]
      data.table::setnames(all_exclusion_wide, "tmp_col", icol)
    }
  }

  all_exclusion_wide[
    j = Status := factor(
      x = data.table::fifelse(
        test = (Sample_Call_Rate + Heterozygosity_Check) > 0,
        yes = "Exclude",
        no = "Check"
      ),
      levels = c("Exclude", "Check")
    )
  ]
  data.table::setcolorder(
    x = all_exclusion_wide,
    neworder = c("FID", "IID", "Status", "Sample_Call_Rate", "Heterozygosity_Check")
  )

  data.table::setorderv(
    x = all_exclusion_wide,
    cols = setdiff(colnames(all_exclusion_wide), c("FID", "IID")),
    order = sign(setdiff(colnames(all_exclusion_wide), c("FID", "IID")) %in% "Status" - 0.5)
  )
}

#' compute_bfile_good_samples
#' @import data.table
compute_bed_good_samples <- function(bfile, exclude, temp_directory, plink) {
  temp_file <- sprintf("%s/data_sample_qc", temp_directory)

  data.table::fwrite(
    x = exclude,
    file = sprintf("%s.txt", temp_file),
    sep = "\t",
    col.names = TRUE,
    quote = FALSE,
    row.names = FALSE
  )

  system(paste(plink,
    "--bfile", bfile,
    "--remove ", sprintf("%s.txt", temp_file),
    "--make-bed",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  out_files <- list.files(
    path = dirname(temp_file),
    pattern = paste(paste0(basename(temp_file), c(".bed", ".bim", ".fam")), collapse = "|"),
    full.names = TRUE
  )

  unlink(
    x = setdiff(
      list.files(
        path = dirname(temp_file),
        pattern = paste0(basename(temp_file), "\\."),
        full.names = TRUE
      ),
      out_files
    ),
    force = TRUE
  )

  out_files
}

#' compute_callrate_snp
#' @import data.table
compute_callrate_snp <- function(bfile, callrate, plink) {
  temp_file <- tempfile(pattern = "crsnp")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))

  system(paste(plink,
    "--bfile", bfile,
    "--missing",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  data.table::fread(file = sprintf("%s.lmiss", temp_file))[
    j = labs := data.table::fifelse(is.na(F_MISS) | F_MISS > 1 - callrate, SNP, NA_character_)
  ][
    F_MISS >= 0.001
  ][
    order(F_MISS)
  ][
    j = status := data.table::fifelse(is.na(labs), "clear", "exclude")
  ]
}

#' plot_callrate_snp
#' @import data.table
#' @import ggplot2
#' @import scales
plot_callrate_snp <- function(data, callrate, max_labels) {
  plot_callrate(
    data = data,
    callrate = callrate,
    max_labels = max_labels
  ) +
    ggplot2::scale_x_continuous(
      labels = scales::comma_format(accurracy = 1),
      trans = "log10",
      guide = ggplot2::guide_axis(check.overlap = TRUE),
      sec.axis = ggplot2::dup_axis(
        labels = function(x) {
          paste0("<b style = 'color:#b22222;'>", scales::comma_format(accurracy = 1)(x), "</b>")
        },
        breaks = data[status %in% "exclude", .N],
        name = NULL
      )
    )
}

#' compute_duplicated_snp
#' @import data.table
compute_duplicated_snp <- function(bfile) {
  list_snps <- data.table::fread(
    file = sprintf("%s.bim", bfile),
    select = c(1, 2, 4),
    col.names = c("CHR", "SNP", "POS")
  )[
    j = CHRPOS := paste(CHR, POS, sep = ":")
  ]

  list_snps[CHRPOS %in% CHRPOS[duplicated(CHRPOS)]][j = status := "exclude"]
}

#' compute_hwe_snp
#' @import data.table
compute_hwe_snp <- function(bfile, hwe, plink) {
  temp_file <- tempfile(pattern = "hwesnp")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))

  system(paste(plink,
    "--bfile", bfile,
    "--hardy",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  data.table::fread(file = sprintf("%s.hwe", temp_file))[
    j = labs := data.table::fifelse(is.na(P) | P < hwe, SNP, NA_character_)
  ][
    order(P)
  ][
    j = status := data.table::fifelse(is.na(labs), "clear", "exclude")
  ]
}

#' plot_hwe_snp
#' @import data.table
#' @import ggplot2
#' @import scales
plot_hwe_snp <- function(data, hwe) {
  ggplot2::ggplot(data = data[P < hwe * 10]) +
    ggplot2::aes(x = seq_along(P), y = P) +
    ggplot2::geom_point(
      colour = scales::viridis_pal(begin = 0.5, end = 0.5)(1),
      shape = 1,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = hwe, colour = "#b22222", linetype = 2) +
    ggplot2::geom_vline(xintercept = data[P < hwe, .N], colour = "#b22222", linetype = 2) +
    ggplot2::scale_x_continuous(
      labels = scales::comma_format(accurracy = 1),
      trans = "log10",
      guide = ggplot2::guide_axis(check.overlap = TRUE),
      sec.axis = ggplot2::dup_axis(
        labels = function(x) {
          paste0("<b style = 'color:#b22222;'>", scales::comma_format(accurracy = 1)(x), "</b>")
        },
        breaks = data[status %in% "exclude", .N],
        name = NULL
      )
    ) +
    ggplot2::scale_y_continuous(trans = pval_trans(alpha = hwe, md = TRUE)) +
    ggplot2::labs(
      x = "Number of Variants",
      y = "Hardy-Weinberg Equilibrium P-value",
      title = "Distribution of Hardy-Weinberg Equilibrium (HWE) P-values",
      subtitle = paste(
        "The horizontal red line (\"<b style = 'color:#b22222;'>---</b>\") denotes the",
        "cut-off for significance above which variants are excluded."
      ),
      caption = paste(
        "Variants with a HWE p-value ≥ 10<sup>-3</sup> are not shown."
      )
    )
}

#' compute_maf_snp
#' @import data.table
compute_maf_snp <- function(bfile, maf, plink) {
  temp_file <- tempfile(pattern = "frqsnp")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))

  system(paste(plink,
    "--bfile", bfile,
    "--freq",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  data.table::fread(file = sprintf("%s.frq", temp_file))[
    j = labs := data.table::fifelse(is.na(MAF) | MAF < maf, SNP, NA_character_)
  ][
    order(MAF)
  ]
}

#' plot_maf_snp
#' @import data.table
#' @import ggplot2
#' @import scales
plot_maf_snp <- function(data) {
  ggplot2::ggplot(
    data = data[!is.na(MAF)][
      j = maf_bin := cut(
        x = MAF,
        breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
        include.lowest = TRUE
      )
    ]
  ) +
    ggplot2::aes(x = maf_bin) +
    ggplot2::geom_bar(fill = scales::viridis_pal(begin = 0.5, end = 0.5)(1)) +
    ggplot2::scale_y_continuous(
      labels = scales::comma_format(accurracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::labs(
      x = "MAF",
      y = "Number of Variants",
      title = "Distribution of Variants Number per Class of Minor Allele Frequency (MAF)"
    )
}

#' compute_variant_qc
#' @import data.table
compute_bed_good_variants <- function(bfile, exclude, temp_directory, project, plink) {
  temp_file <- tempfile(pattern = "data_sample_variant_qc")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))

  out_files_prefix <- sprintf(
    "%s/%s_qc",
    sub("/tmp$", "/plink_qc", temp_directory),
    tolower(project)
  )

  data.table::fwrite(
    x = unique(data.table::rbindlist(exclude, use.names = TRUE, fill = TRUE)[
      i = status %in% "exclude",
      j = c("CHR", "SNP")
    ]),
    file = sprintf("%s.txt", temp_file),
    sep = "\t",
    col.names = TRUE,
    quote = FALSE,
    row.names = FALSE
  )

  system(paste(plink,
    "--bfile", bfile,
    "--set-hh-missing",
    "--exclude", sprintf("%s.txt", temp_file),
    "--make-bed",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)


  if (file.exists(sprintf("%s.hh", temp_file))) {
    system(paste(plink,
      "--bfile", temp_file,
      "--exclude", sprintf("%s.hh", temp_file),
      "--make-bed",
      "--out", out_files_prefix
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  } else {
    system(paste(plink,
      "--bfile", temp_file,
      "--make-bed",
      "--out", out_files_prefix
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  }

  out_files <- list.files(
    path = dirname(out_files_prefix),
    pattern = paste(paste0(basename(out_files_prefix), c(".bed", ".bim", ".fam")), collapse = "|"),
    full.names = TRUE
  )

  unlink(
    x = setdiff(
      list.files(
        path = dirname(out_files_prefix),
        pattern = paste0(basename(out_files_prefix), "\\."),
        full.names = TRUE
      ),
      out_files
    ),
    force = TRUE
  )

  out_files
}

#' compute_vcf_imputation
#' @import data.table
compute_vcf_imputation <- function(
  bfile,
  ref,
  ref_panel,
  ref1kg_fasta,
  perl_script,
  temp_directory,
  project,
  plink,
  bcftools
) {
  temp_file <- tempfile(pattern = "forimputation")
  on.exit(unlink(
    x = list.files(path = dirname(temp_file), pattern = basename(temp_file), full.names = TRUE),
    force = TRUE
  ))

  system(paste(plink,
    "--bfile", bfile,
    "--make-bed",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  system(paste(plink,
    "--bfile", temp_file,
    "--freq",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)

  system(paste(
    "(cd", dirname(temp_file), "&&",
      "perl", perl_script,
      "-b", sprintf("%s.bim", temp_file),
      "-f", sprintf("%s.frq", temp_file),
      "-r", ref_panel,
      switch(EXPR = ref, "HRC" = "-h", "1KG" = "-g --pop ALL"),
      "-l", plink,
    ")"
  ))
  system(sprintf("cd %s && bash Run-plink.sh", dirname(temp_file)))

  chr_exists <- file.exists(sprintf("%s-updated-chr%d.bim", temp_file, 1:23))
  data.table::fwrite(
    x = data.table::data.table(sprintf("%s-updated-chr%d", temp_file, 1:23)[which(chr_exists)]),
    file = sprintf("%s.txt", temp_file),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  system(paste(plink,
    "--merge-list", sprintf("%s.txt", temp_file),
    if (
      anyDuplicated(data.table::fread(sprintf("%s.fam", bfile), header = FALSE, select = 2)) == 0
    ) {
      "--recode vcf-iid"
    } else {
      "--recode vcf"
    },
    "bgz",
    "--out", temp_file
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(sprintf("%s index %s.vcf.gz", bcftools, temp_file))

  data.table::fwrite(
    x = data.table::data.table(
      V1 = 1L:26L,
      V2 = c(1L:22L, "X", "Y", "X", "MT")
    ),
    file = file.path(dirname(temp_file), "plink2ensembl.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  out_file <- sprintf(
    "%s/%s_qc.vcf.gz",
    sub("/tmp$", "/vcf_qc", temp_directory),
    tolower(project)
  )

  if (is.null(ref1kg_fasta)) {
    warning("Reference alleles were not checked, thus the Sanger Imputation Service can't be use!")
    system(sprintf(
      "%s annotate --rename-chrs %s/plink2ensembl.txt %s.vcf.gz -Oz > %s",
      bcftools,
      dirname(temp_file),
      temp_file,
      out_file
    ))
  } else {
    out_tmp_file <- sub(".vcf.gz", "_tmp.vcf.gz", out_file)
    on.exit(unlink(
      x = list.files(path = dirname(out_tmp_file), pattern = basename(out_tmp_file), full.names = TRUE),
      force = TRUE
    ), add = TRUE)
    system(sprintf(
      "%s annotate --rename-chrs %s/plink2ensembl.txt %s.vcf.gz -Oz > %s",
      bcftools,
      dirname(temp_file),
      temp_file,
      out_tmp_file
    ))
    system(sprintf("%s index %s", bcftools, out_tmp_file))

    system(sprintf(
      "%s +fixref %s -Ob -o %s -- -d -f %s -m flip",
      bcftools,
      out_tmp_file,
      out_file,
      ref1kg_fasta
    ))
  }

  system(sprintf("%s index %s", bcftools, out_file))

  list.files(
    path = dirname(out_file),
    pattern = basename(out_file),
    full.names = TRUE
  )
}

#' compute_vcf_dim
#' @import data.table
compute_vcf_dim <- function(vcf) {
  tmp <- data.table::fread(vcf, skip = "#CHROM")
  c(variants = nrow(tmp), samples = ncol(tmp) - grep("FORMAT", colnames(tmp)))
}

#' create_ga_export_directory
#'
create_ga_export_directory <- function(path, project, array) {
  array_directory <- file.path(path, project, array)
  # unlink(x = array_directory, force = TRUE, recursive = TRUE)
  invisible(
    sapply(
      X = file.path(array_directory, c("plink", "vcf", "vcf_imputed_grch37", "vcf_imputed_grch38")),
      FUN = dir.create,
      showWarnings = FALSE,
      recursive = TRUE,
      mode = "0775"
    )
  )
  array_directory
}

#' list_imputed_vcf
#'
list_imputed_vcf <- function(path) {
  out <- list.files(
    path = file.path(path, "vcf_imputed_grch37"),
    pattern = "\\.vcf.gz$",
    full.names = TRUE
  )

  names(out) <- sub(".pbwt_reference_impute.vcf.gz$", "", basename(out))

  out
}

#' compute_vcf_imputed_qc
#' @import data.table
#' @import future.apply
compute_vcf_imputed_qc <- function(vcf, vcftools, uptodate) {
  stopifnot(uptodate)

  temp_dir <- file.path(tempdir(), "imputedqc")
  dir.create(temp_dir, showWarnings = FALSE, mode = "775")
  on.exit(unlink(x = temp_dir, force = TRUE, recursive = TRUE))

  fct_explicit_na <- function(f) {
    factor(
      x = ifelse(is.na(f), "(Missing)", f),
      levels = c(seq_len(nlevels(f)), "(Missing)"),
      labels = c(levels(f), "(Missing)")
    )
  }

  data.table::rbindlist(future.apply::future_lapply(
    X = vcf,
    FUN = function(ichr) {
      invisible(capture.output({
        system(paste(
          vcftools,
          "--gzvcf", ichr,
          "--get-INFO 'INFO'",
          "--get-INFO 'RefPanelAF'",
          "--get-INFO 'AC'",
          "--get-INFO 'AN'",
          "--temp", temp_dir,
          "--out", file.path(temp_dir, basename(ichr))
        ), intern = TRUE, ignore.stdout = TRUE)
      }))
      on.exit(unlink(file.path(temp_dir, basename(ichr)), force = TRUE))

      data.table::setnames(
        x = data.table::fread(
          file = file.path(temp_dir, paste0(basename(ichr), ".INFO")),
          header = TRUE,
          colClasses = c(
            "CHROM" = "character",
            "POS" = "integer",
            "INFO" = "numeric",
            "RefPanelAF" = "numeric",
            "AC" = "numeric",
            "AN" = "numeric"
          ),
          showProgress = FALSE
        ),
        old = "CHROM",
        new = "CHR"
      )[j = AF := AC / AN][
        j = `:=`(
          "bin_af" = fct_explicit_na(
            f = cut(AF, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
          ),
          "bin_info" = fct_explicit_na(
            f = cut(INFO, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
          ),
          "vcf_chr" = sub("^([^.]+)\\..*", "\\1", basename(ichr))
        )
      ]
    }
  ))
}

#' compute_vcf_imputed_qc_info
#' @import data.table
#' @import scales
compute_vcf_imputed_qc_info <- function(data, uptodate) {
  stopifnot(uptodate)
  out <- data[j = list(n = .N), by = bin_info][j = p := n / sum(n)]
  out <- rbind(out, data.table::data.table(bin_info = "Total", n = sum(out[["n"]]), p = 1))
  out[order(bin_info)]
}

#' compute_vcf_imputed_qc_af
#' @import data.table
#' @import scales
compute_vcf_imputed_qc_af <- function(data, uptodate) {
  stopifnot(uptodate)
  out <- data[j = list(n = .N), by = bin_af][j = p := n / sum(n)]
  out <- rbind(out, data.table::data.table(bin_af = "Total", n = sum(out[["n"]]), p = 1))
  out[order(bin_af)]
}

#' save_plot_vcf_imputed_qc
#' @import data.table
#' @import ggplot2
#' @import patchwork
#' @import utils
#' @import ragg
save_plot_vcf_imputed_qc <- function(data, path, uptodate) {
  stopifnot(uptodate)

  temp_file <- file.path(tempdir(), "vcf_imputed_grch37")
  dir.create(path = temp_file, showWarnings = FALSE)
  on.exit(unlink(x = temp_file, force = TRUE, recursive = TRUE))

  my_theme <- ggplot2::theme_minimal(base_family = "Verdana") +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.subtitle = ggplot2::element_text(face = "italic", size = ggplot2::rel(0.80)),
      plot.caption = ggplot2::element_text(face = "italic", size = ggplot2::rel(0.65))
    )

  imp_qc_figure <- function(data, theme) {
    p1 <- ggplot2::ggplot(data = data[!is.na(INFO)]) +
      ggplot2::aes(x = POS, y = seq_along(POS)) +
      ggplot2::geom_point(
        size = 0.1,
        colour = scales::viridis_pal(begin = 0.5, end = 0.5)(1),
        na.rm = TRUE
      ) +
      ggplot2::scale_x_continuous(labels = scales::comma_format(scale = 1 / 1e6)) +
      ggplot2::scale_y_continuous(labels = scales::comma_format(scale = 1 / 1e3)) +
      ggplot2::labs(x = "Position (Mb)", y = "Line Number (x 1,000)")

    p2 <- ggplot2::ggplot(data = data) +
      ggplot2::aes(x = RefPanelAF, y = AF, colour = INFO > 0.8) +
      ggplot2::geom_point(size = 0.1, na.rm = TRUE) +
      ggplot2::scale_colour_viridis_d(begin = 0.2, end = 0.8, na.value = "#b22222") +
      ggplot2::scale_x_continuous(
        labels = scales::percent_format(accuracy = 1, suffix = " %"),
        # expand = c(0, 0),
        breaks = c(0, 0.5, 1),
        limits = c(0, 1)
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1, suffix = " %"),
        # expand = c(0, 0),
        breaks = c(0, 0.5, 1),
        limits = c(0, 1)
      ) +
      ggplot2::labs(x = "HRC Alternate Allele Frequency", y = "Alternate Allele Frequency") +
      ggplot2::guides(colour = "none") # ggplot2::guide_legend(override.aes = list(size = 3)))

    p3 <- ggplot2::ggplot(data = data[!is.na(INFO)]) +
      ggplot2::aes(x = POS, y = INFO, colour = INFO > 0.8) +
      ggplot2::geom_point(size = 0.1, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 0.8, colour = "#b22222", linetype = 2) +
      ggplot2::scale_colour_viridis_d(begin = 0.2, end = 0.8) +
      ggplot2::scale_x_continuous(labels = scales::comma_format(scale = 1 / 1e6)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
      ggplot2::labs(x = "Position (Mb)", y = "INFO Score") +
      ggplot2::guides(colour = "none") # ggplot2::guide_legend(override.aes = list(size = 3)))

    p4 <- ggplot2::ggplot(
      data = data[
        j = list(n = .N),
        by = bin_af
      ][
        j = p := scales::percent(n / sum(n), accuracy = 0.01, suffix = " %")
      ]
    ) +
      ggplot2::aes(x = bin_af, y = n, label = p) +
      ggplot2::geom_bar(fill = scales::viridis_pal(begin = 0.5, end = 0.5)(1), stat = "identity") +
      ggplot2::geom_text(vjust = -0.1, size = 2.5) +
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45))  +
      ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.20))) +
      ggplot2::labs(x = "Alternate Allele Frequency", y = "SNP Count")

    p5 <- ggplot2::ggplot(
      data = data[
        j = list(n = .N),
        by = bin_info
      ][
        j = p := scales::percent(n / sum(n), accuracy = 0.01, suffix = " %")
      ]
    ) +
      ggplot2::aes(x = bin_info, y = n, label = p) +
      ggplot2::geom_bar(fill = scales::viridis_pal(begin = 0.5, end = 0.5)(1), stat = "identity") +
      ggplot2::geom_text(vjust = -0.1, size = 2.5) +
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
      ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.20))) +
      ggplot2::labs(x = "INFO Score", y = "SNP Count")

    lapply(list(p1, p2, p3, p4, p5), `+`, theme)
  }

  file_path <- file.path(path, "vcf_imputed_grch37.zip")

  raster_figures <- sapply(X = c(1:22, "X"), function(ichr) {
    file <- if (ichr == "X") {
      sprintf("%s/chr%s.png", temp_file, ichr)
    } else {
      sprintf("%s/chr%02d.png", temp_file, as.numeric(ichr))
    }

    ragg::agg_png(
      filename = file,
      width = 16, height = 24, units = "cm", res = 300, scaling = 0.75
    )
    print(
      patchwork::wrap_plots(
        imp_qc_figure(data[CHR %in% ichr], my_theme),
        design = "12\n33\n44\n55"
      ) +
        patchwork::plot_annotation(
          title = sprintf("Chromosome %s", ichr),
          tag_levels = "A",
          theme = my_theme
        )
    )
    invisible(dev.off())

    file
  })

  local({
    owd <- getwd()
    setwd(unique(dirname(raster_figures)))
    unlink(file_path)
    utils::zip(zipfile = file_path, files = basename(raster_figures))
    setwd(owd)
  })

  file_path
}

#' is_vcf_imputed_uptodate
#' @import targets
is_vcf_imputed_uptodate <- function(pre, post) {
  max(tar_timestamp_raw(deparse1(substitute(pre)))) <
    min(tar_timestamp_raw(deparse1(substitute(post))))
}

#' compute_related_samples_tab
#' @import data.table
compute_related_samples_tab <- function(relatedness, callrate_samples) {
  samples_bad_callrate <- callrate_samples[!is.na(labs), unique(paste(FID, IID, sep = "_"))]
  merge(
    x = merge(
      x = relatedness[
        i = !paste(FID1, IID1, sep = "_") %in% samples_bad_callrate &
          !paste(FID2, IID2, sep = "_") %in% samples_bad_callrate,
        j = list(FID1, IID1, FID2, IID2, PI_HAT)
      ],
      y = callrate_samples[, list(FID, IID, F_MISS)],
      by.x = c("FID1", "IID1"),
      by.y = c("FID", "IID")
    ),
    y = callrate_samples[, list(FID, IID, F_MISS)],
    by.x = c("FID2", "IID2"),
    by.y = c("FID", "IID"),
    suffixes = c("1", "2")
  )[order(PI_HAT, decreasing = TRUE)]
}

#' save_ga_qc_data
#' @import data.table
save_ga_qc_data <- function(from, report, imputation, exclusion_check, relatedness, ethnicity, to) {
  data.table::fwrite(
    x = exclusion_check,
    file = file.path(to, "quality-control-exclusion-checks.csv")
  )
  data.table::fwrite(
    x = relatedness,
    file = file.path(to, "quality-control-relatedness.csv")
  )
  data.table::fwrite(
    x = ethnicity,
    file = file.path(to, "quality-control-ethnicity.csv")
  )

  all(c(
    file.copy(
      from = list.files(file.path(from, "plink_qc"), full.names = TRUE),
      to = file.path(to, "plink"),
      overwrite = TRUE,
      copy.date = TRUE
    ),
    file.copy(
      from = list.files(file.path(from, "vcf_qc"), full.names = TRUE),
      to = file.path(to, "vcf"),
      overwrite = TRUE,
      copy.date = TRUE
    ),
    file.copy(
      from = imputation,
      to = to,
      overwrite = TRUE,
      copy.date = TRUE
    ),
    file.copy(
      from = report,
      to = file.path(to, "quality-control-report.html"),
      overwrite = TRUE,
      copy.date = TRUE
    ),
    file.exists(file.path(to, "quality-control-exclusion-checks.csv")),
    file.exists(file.path(to, "quality-control-relatedness.csv")),
    file.exists(file.path(to, "quality-control-ethnicity.csv"))
  ))
}
