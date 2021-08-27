#' get_symbol_vep
#' @import data.table
get_symbol_vep <- function(
  input = "snps_locations.txt.gz",
  output_directory = here::here("outputs"),
  genome_assembly = "GRCh38",
  ensembl_version = "104",
  ensembl_species = "homo_sapiens",
  vep_cache = c(
    "server" = "/media/Data/ExternalData/vep_data",
    "docker" = "/disks/DATA/ExternalData/vep_data"
  )
) {
  data.table::fwrite(
    x = unique(data.table::rbindlist(lapply(
      X = list.files(path = input_directory, pattern = "\\.csv.gz$", full.names = TRUE),
      FUN = data.table::fread, select = c("CHR", "BP", "reference_allele", "other_allele")
    )))[
      i = order(CHR, BP),
      j = list(CHR, start = BP, end = BP, alleles = paste(reference_allele, other_allele, sep = "/"), strand = "+")
    ],
    file = file.path(output_directory, "snps_locations.txt.gz"),
    col.names = FALSE, row.names = FALSE, sep = "\t"
  )

  input_docker <- sub(
    "/disks/PROJECT",
    "/media/Datatmp",
    input
  )
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
    '#!/bin/bash',
    '\n\nchmod 777', dirname(input_docker),
    '\n\ndocker run',
    '--rm',
    '--name vep',
    '--volume', paste0(dirname(input_docker), ':/data_dir'),
    paste0('--volume ', vep_cache[["server"]], ':/opt/vep/.vep'),
    paste0('ensemblorg/ensembl-vep:release_', ensembl_version, '.0'),
    '/bin/bash -c "./vep',
    '--input_file', file.path("/data_dir", basename(input_docker)),
    '--cache',
    '--offline',
    '--fork 70',
    '--force_overwrite',
    '--assembly', genome_assembly,
    '--check_existing',
    '--no_check_alleles',
    '--symbol',
    '--output_file', file.path("/data_dir", output_docker),
    '&& cut -f 1-4,13-14', file.path("/data_dir", output_docker),
    '| bgzip --thread 70 -f >', file.path("/data_dir", paste0(output_docker, ".gz")),
    '"',
    '\n\nchmod 775', dirname(input_docker),
    '\n'
  ), file = file.path(output_directory, "run_docker_vep.sh"))

  file.path(output_directory, output_docker)
}

#' format_symbol_vep
#' @import data.table
format_symbol_vep <- function(file) {
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
  vep_annotation <- data.table::fread(file = file,  skip = "#U")[
    j = c("CHR", "POS") := data.table::tstrsplit(Location, ":", fixed = TRUE)
  ][
    j = (c("Gene", "Symbol", "rsid")) :=
      list(
        paste(unique(Gene), collapse = ";"),
        data.table::fifelse(
          test = grepl("SYMBOL=", Extra),
          yes = paste(unique(gsub("^.*SYMBOL=([^;]*);.*$", "\\1", Extra)), collapse = ";"),
          no = NA_character_
        ),
        paste(unique(Existing_variation), collapse = ";")
      ),
    by = "#Uploaded_variation"
  ][j = list(CHR, POS, `#Uploaded_variation`, Gene, Symbol, rsid)]

   data.table::fwrite(
    x = data.table::setnames(unique(vep_annotation), "#Uploaded_variation", "chr_pos_ref_alt"),
    file = sub(".txt.gz", "_formated.txt.gz", file)
  )

  sub(".txt.gz", "_formated.txt.gz", file)
}