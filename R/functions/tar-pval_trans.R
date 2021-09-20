#' pval_trans
#' @import scales
pval_trans <- function(alpha = NULL, md = FALSE, prefix = FALSE, colour = "#b22222") {
  scales::trans_new(
    name = "pval",
    domain = c(0, 1),
    transform = function(x) {
      x[x < .Machine$double.xmin] <- .Machine$double.xmin
      -log(x, 10)
    },
    inverse = function(x) {10^-x},
    breaks = (function(n = 5, digits = 3) {
      function(x) {
        values <- as.numeric(format(c(x, alpha), scientific = TRUE, digits = digits))
        max <- floor(-log(min(values, na.rm = TRUE), base = 10))
        if (max == 0) 1 else sort(unique(c(10^-seq(0, max, by = floor(max / n) + 1), alpha)))
      }
    })(),
    format = (function(x, digits = 3) {
      if (md & nchar(system.file(package = "ggtext")) != 0) {
        x_fmt <- sub(
          "^(.*)e[+]*([-]*)0*(.*)$",
          "\\1 &times; 10<sup>\\2\\3</sup>",
          format(x, scientific = TRUE, digits = digits)
        )
        x_fmt[x %in% c(0, 1)] <- x[x %in% c(0, 1)]
        x_fmt <- sub("^1 &times; ", "", x_fmt)
        if (!is.null(alpha)) {
          alpha_idx <- format(x, scientific = TRUE, digits = digits) ==
            format(alpha, scientific = TRUE, digits = digits)
          x_fmt[alpha_idx] <- paste0(
            "<b style='color:", colour, ";'>",
            if (prefix) "&alpha; = " else "", x_fmt[alpha_idx],
            "</b>"
          )
        }
        x_fmt
      } else {
        x_fmt <- sub(
          "^(.*)e[+]*([-]*)0*(.*)$",
          "\\1 %*% 10^\\2\\3",
          format(x, scientific = TRUE, digits = digits)
        )
        x_fmt[x %in% c(0, 1)] <- x[x %in% c(0, 1)]
        x_fmt <- sub("^1 \\%\\*\\% ", "", x_fmt)
        if (!is.null(alpha)) {
          alpha_idx <- format(x, scientific = TRUE, digits = digits) ==
            format(alpha, scientific = TRUE, digits = digits)
          x_fmt[alpha_idx] <- paste0(if (prefix) "alpha == " else "", x_fmt[alpha_idx])
        }
        parse(text = x_fmt)
      }
    })
  )
}
