#' fortify.manhattan
#' @import data.table
fortify.manhattan <- function(data, x, y, group) {
  map_chro <- c(seq(22), "X", "Y", "X", "Y")
  names(map_chro) <- c(seq(24), "X", "Y")

  `:=` <- data.table::`:=`

  out <- data.table::setnames(
    x = data.table::as.data.table(data),
    old = c(x, y, group),
    new = c("x_pos", "y_pval", "x_chr")
  )
  out[j = x_chr := as.character(x_chr)]
  out[j = x_chr := map_chro[gsub("^chr", "", x_chr, ignore.case = TRUE)]]
  out[j = x_chr := factor(x_chr, levels = intersect(c(seq(22), "X", "Y"), x_chr))]
  out[j = x_pos := as.double(x_pos)]
  out[order(x_chr, x_pos)]
  out[j = x_pos := scales::rescale(x = x_pos, to = c(-0.4, 0.4)), by = "x_chr"]
  out[j = x_pos := x_pos + as.integer(x_chr)]
  data.table::setnames(out, c("x_pos", "y_pval", "x_chr"), c("x", "y", "group"))
}

#' StatManhattan
#' @import ggplot2
StatManhattan <- ggplot2::ggproto("StatManhattan", ggplot2::Stat,
  required_aes = c("x", "y", "group"),
  setup_data = function(data, params) {
    fortify.manhattan(data, "x", "y", "group")
  },
  compute_layer = function(data, scales, params) {
    data
  }
)

#' draw_manhattan
#' @import ggplot2
#' @import data.table
#' @importFrom scales viridis_pal
draw_manhattan <- function(data, x, y, chr, label_y = "P-value", alpha = 0.05) {
  data <- data.table::as.data.table(data)#[, .SD, .SDcols = c(x, y, chr)]
  data.table::setnames(data, c(x, y, chr), c("pos", "pvalue", "chr"), skip_absent = TRUE)
  if (is.numeric(data[["chr"]])) data[, "chr" := lapply(.SD, as.character), .SDcols = "chr"]
  ggplot2::ggplot(data = data)  +
    ggplot2::aes(x = .data[["pos"]], y = .data[["pvalue"]], colour = .data[["chr"]]) +
    ggplot2::geom_point(stat = "manhattan", size = 0.60, na.rm = TRUE) +
    ggplot2::annotate(
      geom = "rect",
      xmin = -Inf, xmax = Inf, ymin = 1, ymax = alpha,
      fill = "#b22222", alpha = 0.2, colour = NA
    ) +
    ggplot2::geom_hline(yintercept = alpha, linetype = 2, colour = "#b22222") +
    ggplot2::scale_x_continuous(
      breaks = 1:24,
      labels = c(1:22, "X", "Y"),
      expand = ggplot2::expansion(add = 0.25)
    ) +
    ggplot2::scale_y_continuous(
      trans = pval_trans(alpha = NULL, md = TRUE, colour = "#b22222"),
      expand = ggplot2::expansion(mult = c(0, 0.2)),
      limits = c(0.05, NA)
    ) +
    ggplot2::scale_colour_manual(values = rep(scales::viridis_pal(begin = 1/4, end = 3/4)(2), 12)) +
    ggplot2::labs(colour = "Chromosome", x = "Chromosome", y = label_y) +
    ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), legend.position = "none")
}
