#' sub_chunk
#' @import knitr
sub_chunk <- function(code, chunk_name = NULL, fig_width = 11.5, fig_height = 5.75) {
  sub_chunk_txt <- sprintf(
    "\n```{r %s, fig.width = %s, fig.height = %s, echo = FALSE}\n(%s)()\n```\n\n",
    chunk_name, fig_width, fig_height, paste0(deparse(function() code), collapse = "")
  )
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk_txt), quiet = TRUE))
}
