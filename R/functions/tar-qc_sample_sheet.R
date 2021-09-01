#' qc_sample_sheet
#' @import data.table
qc_sample_sheet <- function(phenotype, relatedness, ethnicity, methylation) {
  relatedness <- data.table::fread(file = relatedness)

  bad_duplicated_samples <- relatedness[
    i = sub("_.*", "", IID1) == sub("_.*", "", IID2)
  ][
    j = melt.data.table(.SD, measure.vars = paste0("F_MISS", 1:2)),
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

  dt <- phenotype[j = `:=`(
    "Sample_Name" = IID,
    "vcf_id" = IID
  )][# Add related variable
    j = is_related := vcf_id %in% relatedness[
      i = sub("_.*", "", IID1) != sub("_.*", "", IID2),
      j = unique(unlist(.SD)),
      .SDcols = paste0("IID", 1:2)
    ]
  ][# keep best duplicates from genotyping
    i = bad_duplicated_samples,
    j = vcf_id := best_iid,
    on = "IID"
  ]

  dt <- merge(
    x = dt,
    y = data.table::fread(file = ethnicity)[j = -c("cohort")] # Add genetics PCs,
    by.x = "vcf_id",
    by.y = "iid",
    all.x = TRUE
  )

  methylation_sample_sheet <- data.table::fread(methylation)

  merge(# Include methylation sample sheet
    x = dt,
    y = methylation_sample_sheet[
      j = .SD,
      .SDcols = grep(
        pattern = paste(c("^Sample_", "^Sentrix_"), collapse = "|"),
        x = names(methylation_sample_sheet)
      )
    ],
    by = "Sample_Name",
    all = TRUE
  )[
    i = is.na(Sentrix_ID),
    j = Sample_Name := NA_character_
  ][
    i = is.na(PC01),
    j = vcf_id := NA_character_
  ][
    i = is.na(vcf_id),
    j = is_related := NA
  ]
}
