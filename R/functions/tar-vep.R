#' get_variants
#' @import data.table
#' @importFrom future.apply future_lapply
get_variants <- function(path, output_directory = tempdir(), bin_path = list(bcftools = "/usr/bin/bcftools")) {
  if (length(path) == 1 && dir.exists(path)) {
    vcf_files <- list.files(path, pattern = "(\\.vcf.gz|\\.vcf)$", full.names = TRUE)
  } else {
    vcf_files <- path
  }

  unique_snps <- unique(data.table::rbindlist(future.apply::future_lapply(
    X = vcf_files,
    future.globals = FALSE,
    FUN = function(ivcf) {
      xx <- data.table::fread(
        cmd = sprintf(
          paste(
            "%s view %s --min-af 0.05 --exclude 'INFO/INFO < 0.8' --min-alleles 2 --max-alleles 2 --types snps",
            "%s query --format '%%CHROM\t%%POS\t%%REF\t%%ALT\n' --print-header",
            sep = " | "
          ),
          bin_path[["bcftools"]], ivcf, bin_path[["bcftools"]]
        ),
        col.names = c("CHR", "POS", "REF", "ALT")
      )[
        j = list(CHR, start = POS, end = POS, REF = paste(REF, ALT, sep = "/"), strand = "+")
      ]
    }
  )))

  output_file <- file.path(output_directory, "snps_locations.txt.gz")
  fwrite(
    x = unique_snps,
    file = output_file,
    col.names = FALSE, row.names = FALSE, sep = "\t"
  )

  output_file
}

#' get_symbol_vep
#' @import data.table
get_symbol_vep <- function(
  input = "snps_locations.txt.gz",
  output_directory = here::here("outputs", "vep"),
  genome_assembly = "GRCh38",
  ensembl_version = "104",
  ensembl_species = "homo_sapiens",
  vep_cache = c(
    "server" = "/media/Data/ExternalData/vep_data",
    "docker" = "/disks/DATA/ExternalData/vep_data"
  )
) {
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE, mode = "0777")
  }
  input_docker <- sub("/disks/PROJECT", "/media/Datatmp", input)
  output_docker <- paste0("snps_vep_", ensembl_version, ".0_", genome_assembly, ".txt")

  if (
    !file.exists(file.path(
      vep_cache[["docker"]], "homo_sapiens", paste0(ensembl_version, "_", genome_assembly)
    ))
  ) {
    system(sprintf(
      paste(
        "cd %s",
        "curl -sO ftp://ftp.ensembl.org/pub/release-%s/variation/vep/%s_vep_%s_%s.tar.gz",
        "tar xzf %s_vep_%s_%s.tar.gz",
        "rm %s_vep_%s_%s.tar.gz",
        sep = " && "
      ),
      vep_cache[["docker"]],
      ensembl_version, ensembl_species, ensembl_version, genome_assembly,
      ensembl_species, ensembl_version, genome_assembly,
      ensembl_species, ensembl_version, genome_assembly
    ))
  }

  cat(paste(
    "#!/bin/bash",
    "\n\nchmod 777", dirname(input_docker),
    "\n\ndocker run",
    "--rm",
    "--name vep",
    "--volume", paste0(dirname(input_docker), ":/data_dir"),
    paste0("--volume ", vep_cache[["server"]], ":/opt/vep/.vep"),
    paste0("ensemblorg/ensembl-vep:release_", ensembl_version, ".0"),
    '/bin/bash -c "./vep',
    "--input_file", file.path("/data_dir", basename(input_docker)),
    "--cache",
    "--offline",
    "--fork 70",
    "--force_overwrite",
    "--assembly", genome_assembly,
    "--check_existing",
    "--no_check_alleles",
    "--symbol",
    "--output_file", file.path("/data_dir", output_docker),
    "&& cut -f 1-4,13-14", file.path("/data_dir", output_docker),
    "| bgzip --thread 70 -f >", file.path("/data_dir", paste0(output_docker, ".gz")),
    '"',
    "\n\nchmod 775", dirname(input_docker),
    "\n"
  ), file = file.path(output_directory, "run_docker_vep.sh"))

  file.path(output_directory, output_docker)
}

#' format_symbol_vep
#' @import data.table
format_symbol_vep <- function(file, bin_path = list(
    tabix = "/usr/bin/tabix",
    bgzip = "/usr/bin/bgzip"
)) {
  default_file <- file
  file <- paste0(file, ".gz")
  stopifnot(
    file.exists(file) && file.mtime(file.path(dirname(file), "run_docker_vep.sh")) < file.mtime(file)
  )
  if (file.exists(file)) {
    on.exit(unlink(c(
      file.path(dirname(file), "snps_locations.txt.gz"),
      default_file,
      paste0(default_file, "_summary.html")
    )))
  }
  vep_annotation <- unique(data.table::fread(file = file,  skip = "#U"))[
    j = Symbol := data.table::fifelse(
      test = grepl("SYMBOL=", Extra),
      yes = sub("^.*SYMBOL=([^;=]*).*$", "\\1", Extra),
      no = NA_character_
    )
  ][
    j = lapply(.SD, function(x) {
      unique_x <- na.exclude(setdiff(unique(x), ""))
      switch(EXPR = as.character(length(unique_x)),
        "0" = NA_character_,
        "1" = unique_x,
        paste(unique_x, collapse = ";")
      )
    }),
    by = "Location",
    .SDcols = c("Gene", "Existing_variation", "Symbol")
  ][
    j = c("CHR", "POS") := data.table::tstrsplit(Location, ":", fixed = TRUE)
  ][
    j = list(CHROM = CHR, POS, Gene, Symbol, rsid = Existing_variation)
  ]

  data.table::fwrite(
    x = unique(vep_annotation),
    file = sub(".txt.gz", "_formatted.txt", file)
  )
  system(paste(
    bin_path[["bgzip"]], "--force 1", sub(".txt.gz", "_formatted.txt", file),
    "&&",
    bin_path[["tabix"]], "--sequence 1 --begin 2 --end 2", sub(".txt.gz", "_formatted.txt.gz", file)
  ))

  sub(".txt.gz", "_formatted.txt.gz", file)
}
