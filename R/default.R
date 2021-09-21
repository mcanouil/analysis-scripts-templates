message(timestamp(quiet = TRUE))
### Project Setup ==================================================================================
library(here)
project_name <- sub("(.*)_[^_]*\\.Rproj$", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "99-default")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  # library(ragg)
  # library(ggplot2)
  # library(ggtext)
  # library(patchwork)
  # library(data.table)
  # library(future)
  # library(future.callr)
  # library(future.apply)
})


### project setup ==================================================================================
# plan(future.callr::callr, workers = 40)
# message(sprintf("Number of workers: %d", future::nbrOfWorkers()))


### Tables and Figures Theme =======================================================================
# options(
#   ggplot2.discrete.colour = function(...) scale_colour_viridis_d(..., begin = 0.15, end = 0.85),
#   ggplot2.discrete.fill = function(...) scale_fill_viridis_d(..., begin = 0.15, end = 0.85),
#   ggplot2.continuous.colour = function(...) scale_colour_viridis_c(..., begin = 0.15, end = 0.85),
#   ggplot2.continuous.fill = function(...) scale_fill_viridis_c(..., begin = 0.15, end = 0.85)
# )
# theme_set(theme_minimal(base_family = "Verdana"))
# theme_update(
#   plot.title.position = "plot",
#   plot.caption.position = "plot",
#   plot.title = element_markdown(),
#   plot.subtitle = element_markdown(face = "italic"),
#   plot.caption = element_markdown(face = "italic"),
#   axis.title.x = element_markdown(),
#   axis.text.x = element_markdown(),
#   axis.title.y = element_markdown(),
#   axis.text.y = element_markdown()
# )


### Functions ======================================================================================


### Analysis =======================================================================================


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
message(timestamp(quiet = TRUE))
