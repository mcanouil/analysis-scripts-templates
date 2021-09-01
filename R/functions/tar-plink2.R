#' download_plink2
#' @import utils
download_plink2 <- function(url = "http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip", path) {
  plink_path <- file.path(dirname(path), "plink2")
  utils::download.file(url = url, destfile = sprintf("%s.zip", plink_path))
  on.exit(unlink(sprintf("%s.zip", plink_path)))
  utils::unzip(sprintf("%s.zip", plink_path), exdir = dirname(path))
  unlink(sprintf("%s.zip", plink_path))
  Sys.chmod(plink_path, "0777")
  plink_path
}
