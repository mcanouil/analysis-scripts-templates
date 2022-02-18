#' clean_sample_sheet
#'
#' @import here
#' @import readxl
#' @import data.table
clean_sample_sheet <- function(sample_sheet, design = NULL) {
  if (is.null(design)) {
    tmp <- data.table::fread(sample_sheet, skip = "Sample_Name")
  } else {
    tmp <- merge(
      x = data.table::fread(sample_sheet, skip = "Sample_Name"),
      y = unique(data.table::setDT(readxl::read_excel(design))[
        j = c("id", "DNAID", "cohort")
      ]),
      by.x = "Sample_Name",
      by.y = "id",
      all.x = TRUE
    )
  }

  dir.create(
    path = here::here("outputs", "dna_design"),
    recursive = TRUE,
    showWarnings = FALSE,
    mode = "0775"
  )
  data.table::fwrite(
    x = tmp,
    file = here::here("outputs", "dna_design", "epic_design_sheet.csv")
  )

  here::here("outputs", "dna_design", "epic_design_sheet.csv")
}
