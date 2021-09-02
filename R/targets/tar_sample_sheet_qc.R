tar_sample_sheet_qc <- list({
  tar_target(sample_sheet_qc,
    command = qc_sample_sheet(
      phenotype = harmonised_phenotypes,
      exclusion = file.path(ga_export_directory, "quality-control-exclusion-checks.csv"),
      relatedness = file.path(ga_export_directory, "quality-control-relatedness.csv"),
      ethnicity = file.path(ga_export_directory, "quality-control-ethnicity.csv"),
      methylation = ma_csv
    ),
    packages = "data.table"
  )
})