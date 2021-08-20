#' sub_chunk
#'
#' @import knitr
sub_chunk <- function(code, chunk_name = NULL, fig_width = 11.5, fig_height = 5.75) {
  sub_chunk_txt <- paste0(
    "\n```{r ", chunk_name, ", fig.width = ", fig_width, ", fig.height = ", fig_height, ", echo = FALSE}\n",
    "(", paste0(deparse(function() code), collapse = ""), ")()",
    "\n```\n\n"
  )
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk_txt), quiet = TRUE))
}
