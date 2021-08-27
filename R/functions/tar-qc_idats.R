#' Efficiently import idats files mostly using minfi functions.
#'
#' @param directory A `character`. Location of IDAT files, default is the current working directory.
#' @param csv_file A `character`. Path to the sample sheet (csv files) or
#'     name of the sample sheet in `directory`.
#' @param meth_value_type A `character`. Indicates whether you prefer m-values (`"M"`)
#'     or beta-values (`"B"`). Default is `"B"`.
#' @param array_name A `character`. Choose microarray type, eiyther `"450K"` or `"EPIC"`.
#'     Default is `"EPIC"`.
#' @param annotation_version A `character`. Version of the annotation package that should be used.
#'     Default is `"ilm10b4.hg19"` for the `"EPIC"` array
#' @param n_cores An `integer`. The number of cores to use,
#'     i.e., at most how many child processes will be run simultaneously.
#' @param rgSet A `RGChannelSet` object.
#' @param echo A `logical`. Should messages be displayed?
#'
#' @inheritParams qc_idats
#'
#' @import data.table
#' @import ChAMPdata
#'
#' @return A `list`.
#' @export
#'
#' @import data.table
#' @import ENmix
#' @import minfi
#' @import utils
read_idats <- function(
  directory = getwd(),
  csv_file = "csv$",
  meth_value_type = "B",
  filter_beads = TRUE,
  bead_cutoff = 0.05,
  filter_non_cpg = TRUE,
  filter_snps = TRUE,
  population = NULL,
  filter_multihit = TRUE,
  filter_xy = TRUE,
  detection_pvalues = 0.01,
  filter_callrate = TRUE,
  callrate_samples = 0.99,
  callrate_probes = 1,
  norm_background = "oob",
  norm_dye = "RELIC",
  norm_quantile = "quantile1",
  array_name = c("EPIC", "450k"),
  annotation_version = c("ilm10b5.hg38", "ilm10b4.hg19", "ilmn12.hg19"),
  n_cores = 1,
  rgSet = NULL,
  echo = FALSE
) {

  old_r_libs <- Sys.getenv("R_LIBS")
  r_libs <- unique(c(
    .libPaths()[1],
    unlist(strsplit(Sys.getenv("R_LIBS"), .Platform$path.sep))
  ))
  Sys.setenv(R_LIBS = paste(r_libs, collapse = .Platform$path.sep))
  on.exit(Sys.setenv(R_LIBS = old_r_libs))

  if (echo) {
    message(
      "======================\n",
      "Reading IDAT files ...\n",
      "======================",
      appendLF = TRUE
    )
  }

  array_name <- array_name[1]
  annotation_version <- annotation_version[1]

  stopifnot(suppressPackageStartupMessages(
    nchar(system.file(package = "minfi")) > 0 &
      switch(
        EXPR = array_name,
        "450k" = nchar(system.file(package = "IlluminaHumanMethylation450kmanifest")) > 0,
        "EPIC" = nchar(system.file(package = "IlluminaHumanMethylationEPICmanifest")) > 0
      )
  ))

  if (is.null(rgSet) | !inherits(rgSet, "RGChannelSet")) {
    sample_sheet <- read_sample_sheet(directory = directory, csv_file = csv_file, echo = echo)
    rgSet <- read_metharray_exp(sample_sheet = sample_sheet, n_cores = n_cores)
  }
  rgSet@annotation <- switch(
    EXPR = array_name,
    "450k" = c(array = "IlluminaHumanMethylation450k", annotation = annotation_version),
    "EPIC" = c(array = "IlluminaHumanMethylationEPIC", annotation = annotation_version)
  )

  data_detP <- minfi::detectionP(rgSet)
  data_detP[is.na(data_detP)] <- 1

  if (filter_callrate) {
    good_detection <- data_detP < detection_pvalues

    call_rate_samples <- colSums(good_detection) / nrow(good_detection)
    bad_samples <- names(which(call_rate_samples < callrate_samples))

    good_detection <- good_detection[, setdiff(colnames(good_detection), bad_samples)]

    call_rate_cpg <- rowSums(good_detection) / ncol(good_detection)
    bad_cpgs <- names(which(call_rate_cpg < callrate_probes))
  } else {
    bad_samples <- NULL
    bad_cpgs <- NULL
  }

  mset_raw <- mset <- minfi::preprocessRaw(rgSet)

  if (filter_beads) {
    bc <- get_beadcount(rgSet)
    bc2 <- bc[rowSums(is.na(bc)) < bead_cutoff * (ncol(bc)), ]
    mset_f2 <- mset[minfi::featureNames(mset) %in% rownames(bc2), ]
    n_beads_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  if (filter_non_cpg) {
    mset_f2 <- minfi::dropMethylationLoci(mset, dropCH = TRUE)
    n_non_cpg_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  if (filter_snps) {
    ref_population <- c(
      "AFR", "EAS", "EUR",
      "SAS", "AMR", "GWD", "YRI", "TSI", "IBS",
      "CHS", "PUR", "JPT", "GIH", "CHB", "STU",
      "ITU", "LWK", "KHV", "FIN", "ESN", "CEU",
      "PJL", "ACB", "CLM", "CDX", "GBR", "BEB",
      "PEL", "MSL", "MXL", "ASW"
    )

    if (is.null(population) || !(population %in% ref_population)) {
      manifest_hg19 <- switch(
        EXPR = array_name,
        "450k" = get(utils::data("hm450.manifest.hg19", package = "ChAMPdata")),
        "EPIC" = get(utils::data("EPIC.manifest.hg19", package = "ChAMPdata"))
      )
      which_population <- which(manifest_hg19$MASK_general)
    } else {
      manifest_hg19 <- switch(
        EXPR = array_name,
        "450k" = get(utils::data("hm450.manifest.pop.hg19", package = "ChAMPdata")),
        "EPIC" = get(utils::data("EPIC.manifest.pop.hg19", package = "ChAMPdata"))
      )
      which_population <- which(manifest_hg19[, paste("MASK_general", population, sep = "_")])
    }
    maskname <- rownames(manifest_hg19)[which_population]
    mset_f2 <- mset[!minfi::featureNames(mset) %in% maskname, ]
    n_snps_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  if (filter_multihit) {
    multi_hit <- get(utils::data("multi.hit", package = "ChAMPdata"))
    mset_f2 <- mset[!minfi::featureNames(mset) %in% multi_hit$TargetID, ]
    n_multihit_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  if (filter_xy) {
    switch(
      EXPR = array_name,
      "450k" = utils::data("probe.features", package = "ChAMPdata"),
      "EPIC" = utils::data("probe.features.epic", package = "ChAMPdata")
    )
    probe_features <- get("probe.features")
    autosomes <- probe_features[!probe_features$CHR %in% c("X", "Y"), ]
    mset_f2 <- mset[minfi::featureNames(mset) %in% rownames(autosomes), ]
    n_xy_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  mset <- ENmix::preprocessENmix(
    rgSet = rgSet,
    bgParaEst = norm_background,
    dyeCorr = norm_dye,
    QCinfo = NULL,
    exQCsample = FALSE,
    exQCcpg = FALSE,
    exSample = bad_samples,
    exCpG = unique(c(bad_cpgs, setdiff(minfi::featureNames(mset_raw), minfi::featureNames(mset)))),
    nCores = n_cores
  )
  mset <- ENmix::norm.quantile(mdat = mset, method = norm_quantile)

  methylation_matrix <- switch(meth_value_type,
    "B" = minfi::getBeta(mset, "Illumina"),
    "M" = minfi::getM(mset),
    stop('Methylation value type not defined. Only "B" or "M" are available.')
  )

  if (min(methylation_matrix, na.rm = TRUE) <= 0) {
    methylation_matrix[methylation_matrix <= 0] <- min(methylation_matrix[methylation_matrix > 0], na.rm = TRUE)
  }
  if (max(methylation_matrix, na.rm = TRUE) >= 1) {
    methylation_matrix[methylation_matrix >= 1] <- max(methylation_matrix[methylation_matrix < 1], na.rm = TRUE)
  }

  colnames(methylation_matrix) <- minfi::pData(mset)[["Sample_ID"]]
  mset@metadata[[meth_value_type]] <- methylation_matrix
  tmp_phenotypes <- data.table::as.data.table(minfi::pData(rgSet))
  tmp_phenotypes[, "Sample_ID" := lapply(.SD, as.character), .SDcols = "Sample_ID"]
  tmp_means <- colMeans(data_detP)[tmp_phenotypes[["Sample_ID"]]]
  tmp_phenotypes[, "mean_detection_pvalue" := tmp_means]
  tmp_callrate <- (colSums(data_detP < detection_pvalues) / nrow(data_detP))[tmp_phenotypes[["Sample_ID"]]]
  tmp_phenotypes[, "call_rate" := tmp_callrate]
  mset@metadata[["phenotypes"]] <- tmp_phenotypes

  log_msg <- character(0)
  if (filter_callrate) {
    log_msg <- c(log_msg, paste0(
      "Filtering samples with call rate below ",
      paste(format(callrate_samples * 100, digits = 1, nsmall = 1), "%"), ":\n",
      "  - ", format(length(bad_samples), big.mark = ",", digits = 0), " samples were discarded."
    ))
    log_msg <- c(log_msg, paste0(
      "Filtering probes with call rate below ",
      paste(format(callrate_probes * 100, digits = 1, nsmall = 1), "%"), ":\n",
      "  - ", format(length(bad_cpgs), big.mark = ",", digits = 0), " probes were discarded."
    ))
  }
  if (filter_beads) {
    log_msg <- c(log_msg, paste0(
      "Filtering probes with a beadcount lower than three in at least ",
      paste(format(bead_cutoff * 100, digits = 1, nsmall = 1), "%"), " of samples:\n",
      "  - ", n_beads_discarded, " probes were discarded."
    ))
  }
  if (filter_non_cpg) {
    log_msg <- c(log_msg, paste0(
      "Filtering non-cg probes:\n",
      "  - ", n_non_cpg_discarded, " probes were discarded."
    ))
  }
  if (filter_snps) {
    log_msg <- c(log_msg, paste0(
      "Filtering probes with SNPs (Zhou et al., 2016; doi:10.1093/nar/gkw967):\n",
      "  - ", n_snps_discarded, " probes were discarded."
    ))
  }
  if (filter_multihit) {
    log_msg <- c(log_msg, paste0(
      "Filtering probes that align to multiple locations ",
      "(Nordlund et al., 2013; doi:10.1186/gb-2013-14-9-r105):\n",
      "  - ", n_multihit_discarded, " probes were discarded."
    ))
  }
  if (filter_xy) {
    log_msg <- c(log_msg, paste0(
      "Filtering probes on the X or Y chromosome:\n",
      "  - ", n_xy_discarded, " probes were discarded."
    ))
  }
  log_msg <- c(log_msg,
    "Zeros have been replaced with smallest value over zero.",
    "Ones have been replaced with largest value below one.",
    paste0(
      "Data contains:\n",
      "  - ", format(dim(methylation_matrix)[1], big.mark = ",", digits = 0), " probes.\n",
      "  - ", format(dim(methylation_matrix)[2], big.mark = ",", digits = 0), " samples.\n",
      "  - ", format(sum(is.na(methylation_matrix)), big.mark = ",", digits = 0), " missing values."
    )
  )
  if (echo) {
    message(paste(log_msg, collapse = "\n"), appendLF = TRUE)
  }

  list(mset = mset, rgset = rgSet, log = log_msg)
}

#' read_sample_sheet
#'
#' @import data.table
read_sample_sheet <- function(
  directory = ".",
  csv_file = "csv$",
  ignore.case = TRUE,
  recursive = TRUE,
  full.names = TRUE,
  echo = FALSE
) {
  if (file.exists(suppressWarnings(normalizePath(csv_file)))) {
    list_files <- normalizePath(csv_file)
  } else {
    list_files <- list.files(
      path = directory,
      pattern = csv_file,
      full.names = full.names,
      ignore.case = ignore.case,
      recursive = recursive
    )
    if (length(list_files) > 1) {
      warning("More than one CSV file have been found!")
      list_files <- list.files[1]
      if (echo) message("File '", list.files, "' will be used.", appendLF = TRUE)
    }
  }

  data_header <- grep("^\\[DATA\\]", readLines(list_files), ignore.case = TRUE)
  if (length(data_header) == 0) data_header <- 0
  col_names <- colnames(data.table::fread(file = list_files, skip = data_header, nrows = 1))
  default_cols <- c("Sample_ID", "Sentrix_ID", "Sentrix_Position")

  cols_missing <- default_cols[!default_cols %in% col_names]
  if (length(cols_missing) != 0) {
    stop(
      "Sample Sheet must contains the following missing columns:\n",
      "  - ", paste(cols_missing, collapse = "\n  - ")
    )
  }

  sample_sheet <- data.table::fread(file = list_files, skip = data_header)
  data.table::setnames(
    x = sample_sheet,
    old = c("Sample_ID", "Slide", "Array", "Sample_Plate", "Sample_Well"),
    new = c("Sample_ID", "Sentrix_ID", "Sentrix_Position", "Sample_Plate", "Sample_Well"),
    skip_absent = TRUE
  )

  basenames <- sub("_Grn\\.idat.*", "", sapply(
    X = paste0(sample_sheet[["Sentrix_ID"]], "_", sample_sheet[["Sentrix_Position"]], "_Grn.idat"),
    FUN = grep,
    x = list.files(path = directory, recursive = recursive, full.names = TRUE),
    value = TRUE,
    USE.NAMES = FALSE
  ), ignore.case = TRUE)
  sample_sheet[, "Basename" := basenames]

  sample_sheet
}

#' read_metharray
#'
#' @import illuminaio
#' @import minfi
#' @import parallel
read_metharray <- function(files, n_cores) {
  Grn <- Red <- NULL # to avoid note from global variable check
  basenames <- unique(sub("_Grn\\.idat.*|_Red\\.idat.*", "", files))
  for (ichannel in c("Grn", "Red")) {
    i_files <- paste0(basenames, "_", ichannel, ".idat")
    names(i_files) <- basename(basenames)
    i_files_exists <- file.exists(i_files)
    if (!all(i_files_exists)) {
      i_filesgz_exists <- file.exists(paste0(i_files, ".gz"))
      if (!all(i_filesgz_exists)) {
        stop(
           "The following specified files do not exist:\n",
          "  - ", paste(i_files[!i_files_exists], collapse = "\n  - ")
        )
      }
      i_files <- paste0(i_files, ".gz")
    }
    suppressWarnings({
      i_idats <- parallel::mclapply(
        X = i_files, mc.preschedule = FALSE, mc.cores = n_cores,
        FUN = illuminaio::readIDAT
      )
    })
    assign(x = ichannel, value = i_idats)
  }

  G_idats <- lapply(X = Grn, FUN = `[[`, "Quants")
  R_idats <- lapply(X = Red, FUN = `[[`, "Quants")

  common_rows <- as.character(Reduce("intersect", lapply(X = G_idats, FUN = rownames)))

  out <- minfi::RGChannelSetExtended(
    Red = do.call("cbind", lapply(X = R_idats, y = common_rows, FUN = function(x, y) x[y, "Mean"])),
    Green = do.call("cbind", lapply(X = G_idats, y = common_rows, FUN = function(x, y) x[y, "Mean"])),
    RedSD = do.call("cbind", lapply(X = R_idats, y = common_rows, FUN = function(x, y) x[y, "SD"])),
    GreenSD = do.call("cbind", lapply(X = G_idats, y = common_rows, FUN = function(x, y) x[y, "SD"])),
    NBeads = do.call("cbind", lapply(X = G_idats, y = common_rows, FUN = function(x, y) x[y, "NBeads"]))
  )

  rownames(out) <- common_rows

  out
}

#' read_metharray_exp
#'
#' @import methods
#' @import minfi
read_metharray_exp <- function(
  directory = NULL,
  sample_sheet = NULL,
  ignore.case = TRUE,
  recursive = TRUE,
  full.names = TRUE,
  n_cores = 1
) {
  if (is.null(sample_sheet)) {
    if (is.null(directory)) directory <- "."
    green_files <- list.files(
      path = directory,
      pattern = "_Grn.idat.*",
      recursive = recursive,
      ignore.case = ignore.case,
      full.names = full.names
    )
    red_files <- list.files(
      path = directory,
      pattern = "_Red.idat.*",
      recursive = recursive,
      ignore.case = ignore.case,
      full.names = full.names
    )

    if (length(green_files) == 0 || length(red_files) == 0) {
      stop("IDAT files must be provided.")
    }

    common_files <- intersect(sub("_Grn.idat.*", "", green_files), sub("_Red.idat.*", "", red_files))

    if (length(common_files) == 0) {
      stop('"Grn" and "Red" idats files must be provided.')
    }

    common_files_green <- paste0(common_files, "_Grn.idat")
    if (!setequal(common_files_green, green_files)) {
      warning(
        "The following files only exists for the green channel:\n",
        "  - ", paste(setdiff(green_files, common_files_green), collapse = "\n  - ")
      )
    }

    common_files_red <- paste0(common_files, "_Red.idat")
    if (!setequal(common_files_red, red_files)) {
       warning(
        "The following files only exists for the red channel:\n",
        "  - ", paste(setdiff(red_files, common_files_red), collapse = "\n  - ")
      )
    }

    rgSet <- read_metharray(common_files, n_cores)
  } else {
    if (all(!grepl("Basename", names(sample_sheet)))) {
      stop('"Basename" must be provided as a column of "sample_sheet".')
    }

    if (is.null(directory)) {
      files <- sample_sheet[["Basename"]]
    } else {
      files <- file.path(directory, basename(sample_sheet[["Basename"]]))
    }

    rgSet <- read_metharray(files, n_cores)

    pD <- methods::as(sample_sheet, "DataFrame")
    pD[["filenames"]] <- files
    rownames(pD) <- colnames(rgSet)
    rgSet@colData <- pD
  }

  minfi::sampleNames(rgSet) <- rgSet[["Sample_ID"]]

  rgSet
}

#' get_beadcount
#'
#' @import minfi
get_beadcount <- function(x) {
  nb <- minfi::getNBeads(x)
  typeI <- minfi::getProbeInfo(x, type = "I")
  typeII <- minfi::getProbeInfo(x, type = "II")
  locus_names <- minfi::getManifestInfo(x, "locusNames")

  bc_temp <- matrix(
    data = NA_real_,
    ncol = ncol(x),
    nrow = length(locus_names),
    dimnames = list(locus_names, minfi::sampleNames(x))
  )
  bc_temp[typeII$Name, ] <- nb[typeII$AddressA, ]
  bc_temp[typeI$Name, ] <- nb[typeI$AddressB, ]
  bc_temp[typeI$Name, ] <- nb[typeI$AddressA, ]
  bc_temp[which(nb[typeI$AddressA, ] < 3 | nb[typeI$AddressB, ] < 3)] <- NA

  data.frame(bc_temp)
}

#' qc_idats
#'
#' @import data.table
qc_idats <- function(params) {
  sample_sheet <- data.table::fread(params[["csv_file"]], skip = "Sample_Name")
  sample_sheet[
    j = Sample_ID := {
      x <- as.character(1:.N)
      data.table::fifelse(x == "1", as.character(Sample_Name), paste(Sample_Name, x, sep = "_"))
    },
    by = "Sample_Name"
  ]

  if (any(grepl("^[0-9]", sample_sheet[["Sample_ID"]]))) {
    sample_sheet[
      j = Sample_ID := paste0("ID_", Sample_ID)
    ]
  }

  setcolorder(sample_sheet, neworder = "Sample_ID")

  if (!is.null(params[["sex_colname"]])) {
    sample_sheet[
      j = qc_observed_sex := c(
        "1" = 1, "2" = 2, "M" = 1, "F" = 2, "0" = 2
      )[as.character(get(params[["sex_colname"]]))]
    ]
    pca_vars <- intersect(colnames(sample_sheet), unique(c(params[["pca_vars"]], "qc_observed_sex")))
  } else {
    pca_vars <- intersect(colnames(sample_sheet), params[["pca_vars"]])
  }

  data.table::fwrite(x = sample_sheet, file = file.path(tempdir(), "sample_sheet.csv"))

  read_idats(
    directory = params[["data_directory"]],
    csv_file = file.path(tempdir(), "sample_sheet.csv"),
    meth_value_type = "B",
    filter_beads = params[["filter_beads"]],
    bead_cutoff = params[["bead_cutoff"]],
    filter_non_cpg = params[["filter_non_cpg"]],
    filter_snps = params[["filter_snps"]],
    population = params[["population"]],
    filter_multihit = params[["filter_multihit"]],
    filter_xy = params[["filter_xy"]],
    detection_pvalues = params[["detection_pvalues"]],
    filter_callrate = params[["filter_callrate"]],
    callrate_samples = params[["callrate_samples"]],
    callrate_probes = params[["callrate_probes"]],
    norm_background = params[["norm_background"]],
    norm_dye = params[["norm_dye"]],
    norm_quantile = params[["norm_quantile"]],
    array_name = sub(".* ", "", params[["array"]]),
    annotation_version = params[["annotation"]],
    n_cores = min(c(nrow(sample_sheet), 25, detectCores()))
  )
}

#' estimate_cell_composition
#'
#' @import FlowSorted.Blood.EPIC
#' @import minfi
#' @import RefFreeEWAS
#' @import stats
#' @import utils
estimate_cell_composition <- function(data_rgset, data_mset, cell_tissue, array, n_cores) {
  switch(
    EXPR = tolower(cell_tissue),
    "blood" = {
      idol_cpgs <- switch(array,
        "450k" = FlowSorted.Blood.EPIC::IDOLOptimizedCpGs450klegacy,
        "EPIC" = FlowSorted.Blood.EPIC::IDOLOptimizedCpGs
      )
      out <- FlowSorted.Blood.EPIC::estimateCellCounts2(
        rgSet = data_rgset,
        compositeCellType = "Blood",
        processMethod = "preprocessNoob",
        probeSelect = "IDOL",
        cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
        referencePlatform = paste0("IlluminaHumanMethylation", array),
        IDOLOptimizedCpGs = idol_cpgs,
        returnAll = FALSE,
        meanPlot = FALSE,
        verbose = FALSE
      )$counts
      colnames(out) <- paste0("CellT_", colnames(out))
      out
    },
    "cordbloodlegacy" = {
      out <- FlowSorted.Blood.EPIC::estimateCellCounts2(
        rgSet = data_rgset,
        compositeCellType = "CordBlood",
        processMethod = "preprocessNoob",
        probeSelect = "auto",
        cellTypes = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC"),
        referencePlatform = "IlluminaHumanMethylation450k",
        IDOLOptimizedCpGs = NULL,
        returnAll = FALSE,
        meanPlot = FALSE,
        verbose = FALSE
      )$counts
      colnames(out) <- paste0("CellT_", colnames(out))
      out
    },
    "cordblood" = {
      if (nchar(system.file(package = "FlowSorted.CordBloodCombined.450k")) > 0) {
        idol_cpgs <- get(utils::data("IDOLOptimizedCpGsCordBlood", package = "FlowSorted.CordBloodCombined.450k"))
      } else {
        idol_cpgs <- NULL
      }
      # FlowSorted.CordBloodCombined.450k <- FlowSorted.CordBloodCombined.450k::FlowSorted.CordBloodCombined.450k()
      out <- FlowSorted.Blood.EPIC::estimateCellCounts2(
        rgSet = data_rgset,
        compositeCellType = "CordBloodCombined",
        processMethod = "preprocessNoob",
        probeSelect = if (is.null(idol_cpgs)) "auto" else "IDOL",
        cellTypes = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC"),
        referencePlatform = paste0("IlluminaHumanMethylation", array),
        # referenceset = "FlowSorted.CordBloodCombined.450k",
        IDOLOptimizedCpGs = idol_cpgs,
        returnAll = FALSE,
        meanPlot = FALSE,
        verbose = FALSE
      )$counts
      colnames(out) <- paste0("CellT_", colnames(out))
      out
    },
    {
      estimate_k_cluster <- function(Rmat, max_k = 25, n_cores = 1) {
        svdRmat <- RefFreeEWAS::svdSafe(Rmat)
        tmp <- do.call("rbind", mclapply(
          X = 0:max_k,
          mc.cores = n_cores,
          mc.preschedule = FALSE,
          mc_Rmat = Rmat,
          mc_svdRmat = svdRmat,
          FUN = function(Ktest, mc_Rmat, mc_svdRmat) {
            N1 <- dim(mc_Rmat)[1]
            N2 <- dim(mc_Rmat)[2]
            if (Ktest == 0) {
              tmpRminLU <- mc_Rmat
            } else {
              tmpRminLU <- mc_Rmat - mc_svdRmat$u[, 1:Ktest] %*%
                (mc_svdRmat$d[1:Ktest] * t(mc_svdRmat$v[, 1:Ktest]))
            }
            tmpSigSq <- rowSums(tmpRminLU * tmpRminLU) / N2

            c(
              K = Ktest,
              AIC = 2 * (N1 + Ktest * (N1 + N2)) +
                N1 * N2 +
                N2 * sum(log(tmpSigSq)),
              BIC = log(N2) * (N1 + Ktest * (N1 + N2)) +
                N1 * N2 +
                N2 * sum(log(tmpSigSq))
            )
        }))

        list(
          icTable = tmp,
          best = tmp[c(AIC = which.min(tmp[, "AIC"]), BIC = which.min(tmp[, "BIC"])), "K"],
          custom_best = tmp[c(
            AIC = which.max(abs(diff(tmp[, "AIC"])[-1])) + 1,
            BIC = which.max(abs(diff(tmp[, "BIC"])[-1])) + 1
          ), "K"]
        )
      }

      beta_matrix <- stats::na.exclude(minfi::getBeta(data_mset))
      max_k <- min(ncol(beta_matrix), 25)
      k_estimated <- min(estimate_k_cluster(
        Rmat = beta_matrix,
        max_k = max_k,
        n_cores = min(n_cores, max_k)
      )$best)
      mu0 <- RefFreeEWAS::RefFreeCellMixInitialize(
        Y = beta_matrix,
        K = k_estimated,
        Y.Distance = NULL,
        Y.Cluster = NULL,
        largeOK = TRUE,
        dist.method = "euclidean"
      )

      RefFreeCellMixObj <- RefFreeEWAS::RefFreeCellMix(
        Y = beta_matrix,
        mu0 = mu0,
        K = NULL,
        iters = 10,
        Yfinal = NULL,
        verbose = FALSE
      )

      out <- RefFreeCellMixObj[["Omega"]]
      colnames(out) <- paste0("CellT_", 1:ncol(out))
      out
    }
  )
}

#' compute_sex_threshold
#'
#' @import minfi
#' @import stats
compute_sex_threshold <- function(data_rgset, sex_threshold) {
  if (is.null(sex_threshold)) {
    sex_predicted <- minfi::getSex(minfi::mapToGenome(data_rgset), cutoff = -2)
    sex_density <- stats::density(sex_predicted$yMed - sex_predicted$xMed, n = 100000)
    min_diff_xy <- which(diff(sign(diff(sex_density$y))) == 2)
    min_diff_xy <- min_diff_xy[which.min(sex_density$y[min_diff_xy])]
    sex_threshold <- round(x = sex_density$x[min_diff_xy], digits = 3)
  }

  sex_threshold
}

#' check_sex
#'
#' @import data.table
#' @import minfi
check_sex <- function(data_rgset, sex_threshold) {
  sex_predicted <- data.table::as.data.table(
    minfi::getSex(minfi::mapToGenome(data_rgset), cutoff = sex_threshold)
  )[j = Sample_ID := minfi::sampleNames(data_rgset)]
  data.table::setnames(
    x = sex_predicted,
    old = c("xMed", "yMed", "predictedSex"),
    new = paste0("qc_", c("xmedian", "ymedian", "predicted_sex"))
  )
  sex_predicted
}

#' compute_phenotypes
#'
#'  @import data.table
compute_phenotypes <- function(data_mset, cell, sex_predicted) {
  raw_phenotypes <- data.table::as.data.table(data_mset@metadata[["phenotypes"]])[
    j = c("Sample_ID", "Sample_Plate", "Sentrix_ID") := lapply(.SD, as.character),
    .SDcols = c("Sample_ID", "Sample_Plate", "Sentrix_ID")
  ]
  if (inherits(cell, c("data.table", "data.frame", "matrix"))) {
    raw_phenotypes <- merge(
      x = raw_phenotypes,
      y = data.table::as.data.table(cell, keep.rownames = "Sample_ID"),
      by = "Sample_ID"
    )
  }

  if (inherits(sex_predicted, c("data.table", "data.frame", "matrix"))) {
    raw_phenotypes <- merge(
      x = raw_phenotypes,
      y = sex_predicted[
        j = c("Sample_ID", "qc_predicted_sex") := list(
          as.character(Sample_ID),
          c("M" = 1, "F" = 2)[qc_predicted_sex]
        )
      ],
      by = "Sample_ID"
    )[j = qc_sex_discrepancy := is.na(qc_observed_sex) | qc_observed_sex != qc_predicted_sex]

    raw_phenotypes[
      j = c("qc_predicted_sex", "qc_observed_sex") :=
        lapply(.SD, factor, levels = c(1, 2), labels = c("Male", "Female")),
      .SDcols = c("qc_predicted_sex", "qc_observed_sex")
    ]
  }

  raw_phenotypes
}

#' normalise_mset
#'
#' @import ENmix
#' @import minfi
#' @import sva
normalise_mset <- function(data_mset, phenotypes) {
  norm_beta <- ENmix::rcp(mdat = data_mset)
  if (length(unique(minfi::pData(data_mset)[["Sentrix_ID"]])) > 1) {
    norm_beta <- sva::ComBat(
      dat = norm_beta,
      batch = factor(minfi::pData(data_mset)[["Sentrix_ID"]])
    )[rownames(data_mset), ]
  }
  colnames(norm_beta) <- minfi::pData(data_mset)[["Sample_ID"]]
  if (min(norm_beta, na.rm = TRUE) <= 0) {
    norm_beta[norm_beta <= 0] <- min(norm_beta[norm_beta > 0])
  }
  if (max(norm_beta, na.rm = TRUE) >= 1) {
    norm_beta[norm_beta >= 1] <- max(norm_beta[norm_beta < 1])
  }
  data_mset@metadata[grep("_values", names(data_mset@metadata))] <- NULL
  data_mset@metadata[["norm_beta_values"]] <- norm_beta
  data_mset@metadata[["phenotypes"]] <- phenotypes

  data_mset
}

#' mset_pca_plot
#'
#' @import data.table
#' @import flashpcaR
#' @import ggplot2
#' @import ggtext
#' @import minfi
#' @import patchwork
#' @import scales
#' @import stats
#' @import utils
mset_pca_plot <- function(data, normalised_mset, pca_vars) {
  phenotypes <- normalised_mset@metadata[["phenotypes"]]
  pca_vars <- intersect(colnames(phenotypes), unique(c(pca_vars, "qc_observed_sex")))

  lapply(
    X = list(
    	"Raw &beta;-values" = `colnames<-`(minfi::getBeta(data$rgset), minfi::pData(data$rgset)[, "Sample_ID"]),
    	"Normalised &beta;-values" = normalised_mset@metadata[["norm_beta_values"]]
    ),
    FUN = function(data_batch) {
      pca_methylation <- data_batch
      pca_methylation <- pca_methylation[rowSums(is.na(pca_methylation)) == 0, ]
      pca_phenotypes <- phenotypes[Sample_ID %in% colnames(pca_methylation)]

      n_comp <- min(10, ncol(pca_methylation))
      fig_n_comp <- min(3, ncol(pca_methylation))

      keep_technical <- names(which(sapply(pca_phenotypes[
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

      pca_res <- flashpcaR::flashpca(X = t(pca_methylation), stand = "sd", ndim = n_comp)

      pca_dfxy <- data.table::as.data.table(pca_res[["vectors"]], keep.rownames = "Sample_ID")
      data.table::setnames(
        x = pca_dfxy,
        old = setdiff(names(pca_dfxy), "Sample_ID"),
        new = sprintf("PC%02d", as.numeric(gsub("V", "", setdiff(names(pca_dfxy), "Sample_ID"))))
      )
      pca_dfxy <- merge(x = pca_dfxy, y = pca_phenotypes, by = "Sample_ID")

      p_inertia <- ggplot2::ggplot(
        data = data.table::data.table(
          y = pca_res[["pve"]],
          x = sprintf("PC%02d", seq_along(pca_res[["pve"]]))
        )[x %in% sprintf("PC%02d", 1:fig_n_comp)]
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
      )[pc %in% sprintf("PC%02d", 1:n_comp)][
        j = {
          m <- stats::model.matrix(
            object = stats::as.formula(
              object = paste0("values ~ ", paste(keep_technical, collapse = " + "))
            ),
            data = .SD
          )

          if (qr(m)$rank == ncol(m)) {
            out <- data.table::as.data.table(
              stats::anova(
                stats::lm(
                  formula = stats::as.formula(
                    object = paste0("values ~ ", paste(keep_technical, collapse = " + "))
                  ),
                  data = .SD
                )
              ),
              keep.rownames = "term"
            )[term != "Residuals"]
          } else {
            out <- data.table::rbindlist(
              lapply(X = keep_technical, .data = .SD, FUN = function(.x, .data) {
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
          out[j = full_rank := qr(m)$rank == ncol(m)]
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
              cols = levels(asso_dt[["pc"]])[1:n_comp],
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
                  X = utils::combn(sprintf("PC%02d", 1:fig_n_comp), 2),
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
  )
}

#' plot_callrate_ma
#'
#' @import data.table
#' @import ggplot2
#' @import scales
#' @import ggrepel
plot_callrate_ma <- function(data, callrate, max_labels) {
  ggplot2::ggplot(
    data = data[order(call_rate), list(Sample_ID, call_rate)][
      j = c("x", "labs") := list(1:.N, ifelse(call_rate < callrate, Sample_ID, NA))
    ]
  ) +
    ggplot2::aes(x = seq_along(call_rate), y = call_rate) +
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
    ggplot2::labs(x = "Number of Samples", y = "Call Rate", title = "Sample Call Rate")
}

#' plot_check_methylation_sex
#'
#' @import data.table
#' @import ggplot2
#' @import scales
#' @import ggrepel
#' @import patchwork
plot_check_methylation_sex <- function(data, sex_threshold) {
  axis_limits <- range(data[j = c("qc_xmedian", "qc_ymedian")], na.rm = TRUE)

  p0 <- ggplot2::ggplot(data = data) +
    ggplot2::aes(x = qc_ymedian - qc_xmedian) +
    ggplot2::geom_density(na.rm = TRUE) +
    ggplot2::geom_vline(
      xintercept = sex_threshold,
      linetype = 2,
      colour = "#b22222",
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      sec.axis = ggplot2::dup_axis(
        name = NULL,
        breaks = function(x) unique(c(scales::breaks_extended()(x), sex_threshold)),
        labels = function(x) ifelse(x == sex_threshold, paste0("<b style='color:#b22222;'>", x, "</b>"), "")
      )
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      x = paste0(
        "Y Chromosome Median Total Intensity (log<sub>2</sub>)<br>",
        "- X Chromosome Median Total Intensity (log<sub>2</sub>)"
      ),
      y = "Density"#,
      # title = "Sex Threshold Detection"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::calc_element(element = "text", ggplot2::theme_get())$size / 3),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA),
      panel.background = ggplot2::element_rect(fill = NA, colour = NA),
      axis.ticks.length = ggplot2::unit(0, "pt"),
      plot.margin = ggplot2::margin(10, 5.5, 5.5, 5.5, unit = "pt")
    )

  p <- ggplot2::ggplot(data = data) +
    ggplot2::aes(
      x = qc_xmedian,
      y = qc_ymedian,
      shape = factor(qc_observed_sex),
      colour = factor(qc_observed_sex)
    ) +
    ggplot2::geom_polygon(
      data = data.frame(
        qc_xmedian = c(c(0, 0, 20), c(0, 20, 20)),
        qc_ymedian = c(c(0, 20, 20), c(0, 0, 20)) + sex_threshold,
        qc_predicted_sex = factor(rep(c(1, 2), each = 3), levels = c(1, 2), labels = c("Male", "Female"))
      ),
      mapping = ggplot2::aes(x = qc_xmedian, y = qc_ymedian, fill = qc_predicted_sex),
      alpha = 0.1,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_abline(
      data = data.frame(Threshold = paste("=", sex_threshold), Seuil = sex_threshold),
      mapping = ggplot2::aes(intercept = Seuil, slope = 1, linetype = Threshold),
      colour = "#b22222",
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = ~ .x[(!qc_sex_discrepancy)],
      size = 2,
      na.rm = TRUE
    ) +
    ggplot2::stat_ellipse(data = ~ .x[(!qc_sex_discrepancy)], na.rm = TRUE, show.legend = FALSE) +
    ggplot2::geom_point(
      data = ~ .x[(qc_sex_discrepancy)],
      colour = "#b22222",
      size = 4,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggrepel::geom_label_repel(
      data = ~ .x[(qc_sex_discrepancy)],
      mapping = ggplot2::aes(x = qc_xmedian, y = qc_ymedian, label = Sample_ID),
      segment.colour = "black",
      colour = "black",
      min.segment.length = ggplot2::unit(0, "lines"),
      size = 2,
      nudge_x = -1,
      inherit.aes = FALSE,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::scale_colour_viridis_d(begin = 0.2, end = 0.8, drop = FALSE) +
    ggplot2::scale_fill_viridis_d(begin = 0.2, end = 0.8, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(22, 21), drop = FALSE) +
    ggplot2::scale_linetype_manual(values = 2) +
    ggplot2::labs(
      x = paste("X Chromosome<br><i>Median Total Intensity (log<sub>2</sub>)</i>"),
      y = paste("Y Chromosome<br><i>Median Total Intensity (log<sub>2</sub>)</i>"),
      colour = "Predicted",
      fill = "Predicted",
      shape = "Observed",
      linetype = "Sex Threshold",
      title = "Sex Check Using Methylation Intensity On X/Y Chromosomes"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1)),
      shape = ggplot2::guide_legend(order = 2, override.aes = list(size = 4)),
      linetype = ggplot2::guide_legend(order = 3),
      colour = "none"
    ) +
    ggplot2::theme(
      legend.key.height = ggplot2::unit(1.2, "lines"),
      legend.key.width = ggplot2::unit(1.2, "lines"),
      legend.spacing.y = ggplot2::unit(5.5, "pt")
    ) +
    ggplot2::coord_cartesian(xlim = axis_limits, ylim = axis_limits)

  p + patchwork::inset_element(p0, 0.80, 0, 1, 0.25, align = "full")
}

#' plot_cell_composition
#'
#' @import data.table
#' @import ggplot2
#' @import scales
#' @import ggrepel
#' @import patchwork
#' @import ggdendro
#' @import stats
plot_cell_composition <- function(data, max_labels) {
  cell_cols <- grep("^CellT_", names(data), value = TRUE)
  dd_row <- stats::as.dendrogram(
    stats::hclust(
      d = stats::dist(data[j = ..cell_cols], method = "euclidean"),
      method = "ward.D2"
    )
  )
  dd_col <- stats::as.dendrogram(
    stats::hclust(
      d = stats::dist(data.table::transpose(data[j = ..cell_cols]), method = "euclidean"),
      method = "ward.D2"
    )
  )
  p_heatmap <- list(
    ggplot2::ggplot(
      data = data.table::melt(
        data[j = .SD, .SDcols = c("Sample_ID", cell_cols)],
        measure.vars = grep("^CellT_", names(data), value = TRUE)
      )[
        j = c("Sample_ID", "variable") :=
          list(
            factor(Sample_ID, levels = data[stats::order.dendrogram(dd_row), Sample_ID]),
            factor(variable, levels = cell_cols[stats::order.dendrogram(dd_col)])
          )
      ]
    ) +
      ggplot2::aes(
        x = variable,
        y = Sample_ID,
        fill = scales::rescale(value, to = c(0, 1))
      ) +
      ggplot2::geom_tile() +
      ggplot2::scale_x_discrete(
        expand = c(0, 0),
        labels = function(x) gsub("CellT_", "", x),
        position = "top"
      ) +
      ggplot2::scale_y_discrete(position = "right", expand = c(0, 0)) +
      ggplot2::scale_fill_viridis_c(
        limits = c(0, 1),
        breaks = c(0, 0.5, 1),
        labels = scales::percent_format(accuracy = 1, suffix = " %"),
        guide = ggplot2::guide_colourbar(
          title = "Composition",
          title.position = "top",
          title.hjust = 0.5,
          direction = "horizontal",
          barwidth = ggplot2::unit(8, units = "lines"),
          raster = TRUE
        )
      ) +
      # theme_minimal() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.y.right = if (nrow(data) > max_labels) ggplot2::element_blank() else ggplot2::element_text(),
        axis.ticks = ggplot2::element_line(colour = "black"),
        axis.ticks.length = ggplot2::unit(x = 0.1, units = "line")
      ) +
      ggplot2::labs(x = "Cell Type", y = "Sample"),

    ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = ggdendro::segment(ggdendro::dendro_data(dd_col, type = "rectangle")),
        mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        size = 0.5
      ) +
      ggplot2::theme_void() +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = c(0.5, 0.5))) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0, 0.1))),

    ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = ggdendro::segment(ggdendro::dendro_data(dd_row, type = "rectangle")),
        mapping = ggplot2::aes(x = y, y = x, xend = yend, yend = xend),
        size = 0.5
      ) +
      ggplot2::theme_void() +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0.5, 0.5))),

    patchwork::guide_area()
  )

  patchwork::wrap_plots(
    p_heatmap, design = "BD\nAC", guides = "collect", widths = c(2/3, 1/3), heights = c(1/3, 2/3)
  ) +
    patchwork::plot_annotation(
      caption = "IDOL optimised CpGs are used when available in the relevant reference set or method."
    )
}

#' create_ma_export_directory
#'
create_ma_export_directory <- function(path, project, array) {
  array_directory <- file.path(path, project, array)
  # unlink(x = array_directory, force = TRUE, recursive = TRUE)
  dir.create(
    path = array_directory,
    showWarnings = FALSE,
    recursive = TRUE,
    mode = "0775"
  )
  array_directory
}

#' export_ma_data
#'
#' @import data.table
#' @import readr
export_ma_data <- function(data_idats, mset, array, output_directory) {
  readr::write_rds(
    x = data_idats,
    file = file.path(output_directory, paste0(array, "_idats.rds"))
  )
   readr::write_rds(
    x = mset,
    file = file.path(output_directory, paste0(array, "_QC_mset.rds"))
  )
  data.table::fwrite(
    x = data.table::as.data.table(mset@metadata[["norm_beta_values"]], keep.rownames = "cpg_id"),
    file = file.path(output_directory, paste0(array, "_QC_betavalues.csv.gz"))
  )
  data.table::fwrite(
    x = data.table::as.data.table(mset@metadata[["phenotypes"]]),
    file = file.path(output_directory, paste0(array, "_QC_phenotypes.csv"))
  )

  file.path(
    output_directory,
    paste0(array, c("_idats.rds", "_QC_mset.rds", "_QC_betavalues.csv.gz", "_QC_phenotypes.csv"))
  )
}

#' save_ma_qc_data
#'
save_ma_qc_data <- function(from, report, to) {
  all(c(
    file.copy(
      from = from,
      to = to,
      overwrite = TRUE,
      copy.date = TRUE
    ),
    file.copy(
      from = report,
      to = file.path(to, "quality-control-report.html"),
      overwrite = TRUE,
      copy.date = TRUE
    )
  ))
}
