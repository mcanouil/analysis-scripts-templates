#' do_crossmap
#' @importFrom future.apply future_lapply
do_crossmap <- function(
  path = NULL,
  output_directory = NULL,
  ref_fasta = "/disks/DATA/ExternalData/1kg/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
  chain_file = "/disks/DATA/ExternalData/1kg/chain_files/hg19ToHg38.over.chain.gz",
  bin_path = list(
    crossmap = "/usr/local/bin/CrossMap.py",
    bcftools = "/usr/bin/bcftools"
  )
) {
  if (!all(sapply(bin_path, file.exists))) {
    stop('"bin_path" must contains valid named path to "crossmap" and "bcftools"!')
  }

  if (length(path) == 1 && dir.exists(path)) {
    vcf_files <- list.files(path, pattern = "(\\.vcf.gz|\\.vcf)$", full.names = TRUE)
  } else {
    vcf_files <- path
  }

  tmpdir <- file.path(tempdir(), "crossmap_tmp")
  dir.create(path = tmpdir, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE))

  tmp_path_converted <- future.apply::future_lapply(
    X = vcf_files,
    ref_fasta = ref_fasta,
    chain_file = chain_file,
    bin_path = bin_path,
    tmpdir = tmpdir,
    future.globals = FALSE,
    FUN = function(vcf, ref_fasta, chain_file, bin_path, tmpdir) {
      output_file <- file.path(tmpdir, basename(vcf))
      system(sprintf(
        paste(
          "%s vcf %s %s %s %s",
          "%s sort %s --output-type z --output %s",
          "%s index --force %s",
          sep = " && "
        ),
        bin_path[["crossmap"]], chain_file, vcf, ref_fasta, output_file,
        bin_path[["bcftools"]], output_file, output_file,
        bin_path[["bcftools"]], output_file
      ))
      if (file.exists(output_file)) {
        output_file
      } else {
        NULL
      }
    }
  )

  if (length(not_converted <- tmp_path_converted[sapply(tmp_path_converted, is.null)]) > 0) {
    warning(sprintf(
      "Crossmap for the following VCF files resulted in an error:\n%s",
      paste(paste("  *", not_converted), collapse = "\n")
    ))
  }

  combined_vcf <- file.path(tmpdir, "all_chromosomes.vcf.gz")
  system(sprintf(
    paste(
      "%s concat %s --allow-overlaps",
      "%s sort --output-type z --output %s",
      sep = "|"
    ),
    bin_path[["bcftools"]], paste(tmp_path_converted[!sapply(tmp_path_converted, is.null)], collapse = " "),
    bin_path[["bcftools"]], combined_vcf
  ))

  path_converted <- future.apply::future_lapply(
    X = vcf_files,
    vcf_file = combined_vcf,
    bin_path = bin_path,
    output_directory = output_directory,
    future.globals = FALSE,
    FUN = function(chr_vcf_file, vcf_file, bin_path, output_directory) {
      output_file <- file.path(output_directory, basename(chr_vcf_file))
      system(sprintf(
        "%s view %s --regions %s --output-type z --output %s",
        bin_path[["bcftools"]],
        paste(paste0(c("", "chr"), sub("\\..*", "", chr_vcf_file)), collapse = ","),
        output_file
      ))
      output_file
    }
  )

  if (length(not_converted <- path_converted[sapply(path_converted, is.null)]) > 0) {
    warning(sprintf(
      "The following final crossmapped VCF files were not produced:\n%s",
      paste(paste("  *", path_converted), collapse = "\n")
    ))
  }

  path_converted
}
