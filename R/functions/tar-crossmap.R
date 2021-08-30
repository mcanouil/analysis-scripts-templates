#' do_crossmap
#' @importFrom future.apply future_lapply future_apply
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
  if (!dir.exists(output_directory)) {
    message(sprintf("Directory %s did not exist and was created!", output_directory))
    dir.create(output_directory, recursive = TRUE)
  }

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

  chromosomes <- unique(sub("\\..*", "", basename(vcf_files)))

  tmp_path_converted <- future.apply::future_lapply(
    X = vcf_files,
    ref_fasta = ref_fasta,
    chain_file = chain_file,
    bin_path = bin_path,
    tmpdir = tmpdir,
    chromosomes = chromosomes,
    future.globals = FALSE,
    FUN = function(vcf, ref_fasta, chain_file, bin_path, tmpdir, chromosomes) {
      temp_file <- file.path(tmpdir, basename(vcf))
      system(sprintf(
        paste(
          "%s vcf %s %s %s %s",
          "%s sort %s --output-type z --output %s",
          "%s index --force %s",
          sep = " && "
        ),
        bin_path[["crossmap"]], chain_file, vcf, ref_fasta, temp_file,
        bin_path[["bcftools"]], temp_file, temp_file,
        bin_path[["bcftools"]], temp_file
      ))

      if (!file.exists(temp_file)) {
        return(NULL)
      }

      output_files <- sprintf("%s/truechr_%s__%s", tmpdir, chromosomes, basename(vcf))
      names(output_files) <- chromosomes

      for (chr in chromosomes) {
        system(sprintf(
          paste(
           "%s view %s --regions %s --output-type z --output %s",
           "%s index --force %s",
            sep = " && "
          ),
          bin_path[["bcftools"]],
            temp_file,
            paste(paste0(c("", "chr"), chr), collapse = ","),
            output_files[[chr]],
          bin_path[["bcftools"]],
            output_files[[chr]]
        ))
      }

      output_files
    }
  )

  if (length(not_converted <- vcf_files[sapply(tmp_path_converted, is.null)]) > 0) {
    warning(sprintf(
      "Crossmap for the following VCF files resulted in an error:\n%s",
      paste(paste("  *", not_converted), collapse = "\n")
    ))
  }

  path_converted <- future.apply::future_apply(
    X = do.call("rbind", tmp_path_converted),
    MARGIN = 2,
    bin_path = bin_path,
    output_directory = output_directory,
    future.globals = FALSE,
    FUN = function(chr_vcf_files, bin_path, output_directory) {
      combined_chr_vcf <- sprintf(
        "%s/%s.vcf.gz",
        output_directory,
        unique(sub("truechr_([^_]+)__.*", "\\1", basename(chr_vcf_files)))
      )
      system(sprintf(
        paste(
          "%s concat %s --allow-overlaps | %s sort --output-type z --output %s",
          "%s index --force %s",
          sep = " && "
        ),
        bin_path[["bcftools"]], paste(chr_vcf_files, collapse = " "),
        bin_path[["bcftools"]], combined_chr_vcf,
        bin_path[["bcftools"]], combined_chr_vcf
      ))

      if (!file.exists(combined_chr_vcf)) {
        return(NULL)
      }

      combined_chr_vcf
    }
  )

  if (length(not_converted <- path_converted[sapply(path_converted, is.null)]) > 0) {
    warning(sprintf(
      "The following final crossmapped VCF files were not produced:\n%s",
      paste(paste("  *", path_converted), collapse = "\n")
    ))
  }

  unname(path_converted)
}
