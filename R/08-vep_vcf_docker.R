### Project Setup ==================================================================================
library(here)
project_name <- gsub("(.*)_.*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "08-vep_vcf_docker")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

dir.create(file.path(output_directory, "vcfs_qc"), recursive = TRUE, showWarnings = FALSE, mode = "0775")

input_directory <- file.path("/disks/DATA/Projects", project_name, "QC")
server_directory <- file.path("/media/Datatmp", project_name)
vep_cache <- c(
  "server" = "/media/Data/ExternalData/vep_data", 
  "docker" = "/disks/DATA/ExternalData/vep_data"
)
genome_assembly <- "GRCh38"
ensembl_version <- "102"
ensembl_species <- "homo_sapiens"


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})


### Compile SNPs List ==============================================================================
if (!file.exists(file.path(output_directory, "snps_locations.txt.gz"))) {
  cat(
    fread(file.path(input_directory, "PreciNASH_omics_samples_qc.csv.gz"))[["omni"]],
    file = file.path(output_directory, "samples_to_keep.txt"),
    sep = "\n"
  )
  
  vcfs <- list.files(
    path = file.path(input_directory, "Omni2.5", "vcf_imputed_hg38_ucsc"),
    pattern = "[^X].vcf.gz$",
    full.names = TRUE
  )
  names(vcfs) <- sprintf("chr%02d", as.numeric(gsub(".pbwt_reference_impute.vcf.gz$", "", basename(vcfs))))
  
  unique_snps <- unique(rbindlist(mclapply(
    X = names(vcfs),
    mc.cores = 11,
    mc.preschedule = FALSE,
    FUN = function(ivcf) {
      fwrite(
        x = fread(
          cmd = paste(
            "vcftools --gzvcf", vcfs[ivcf],
            "--get-INFO 'INFO'",
            "--stdout"
          ),
          header = TRUE,
          colClasses = c("INFO" = "numeric")
        )[INFO < 0.8, c(1, 2)],
        file = file.path(output_directory, "vcfs_qc", paste0(ivcf, "_lowqual.txt")),
        sep = "\t"
      )
  
      system(paste(
        "vcftools --gzvcf", vcfs[ivcf],
        "--keep", file.path(output_directory, "samples_to_keep.txt"),
        "--exclude-positions", file.path(output_directory, "vcfs_qc", paste0(ivcf, "_lowqual.txt")),
        "--remove-indels",
        "--remove-filtered-all",
        "--recode-INFO-all",
        "--maf 0.05",
        "--hwe 0.005",
        "--recode",
        "--out", file.path(output_directory, "vcfs_qc", ivcf)
      ))
      system(paste("bgzip -f", file.path(output_directory, "vcfs_qc", paste0(ivcf, ".recode.vcf"))))
      system(paste("tabix -pvcf -f", file.path(output_directory, "vcfs_qc", paste0(ivcf, ".recode.vcf.gz"))))
  
      system(paste(
        "bcftools annotate",
        "--set-id +'%CHROM:%POS'",
        "--output-type z",
        "--output", file.path(output_directory, "vcfs_qc", paste0(ivcf, "_qc.vcf.gz")),
        file.path(output_directory, "vcfs_qc", paste0(ivcf, ".recode.vcf.gz"))
      ))
      system(paste("tabix -pvcf -f", file.path(output_directory, "vcfs_qc", paste0(ivcf, "_qc.vcf.gz"))))
  
      unlink(c(
        file.path(output_directory, "vcfs_qc", paste0(ivcf, "_lowqual.txt")),
        file.path(output_directory, "vcfs_qc", paste0(ivcf, ".recode.vcf.gz")),
        file.path(output_directory, "vcfs_qc", paste0(ivcf, ".recode.vcf.gz.tbi"))
      ))
        
      fread(
        cmd = paste(
          "vcftools --gzvcf", vcfs[ivcf],
          "--get-INFO 'INFO'",
          "--stdout"
        ),
        header = TRUE,
        colClasses = c("INFO" = "numeric")
      )[INFO >= 0.8][
        j = (c("REF", "strand")) := list(paste(REF, ALT, sep = "/"), "+")
      ][
        j = list(CHR = CHROM, start = POS, end = POS, REF, strand)
      ]
    }
  )))
  
  fwrite(
    x = unique_snps, 
    file = file.path(output_directory, "snps_locations.txt.gz"), 
    col.names = FALSE, row.names = FALSE, sep = "\t"
  )
}


### Docker Command =================================================================================
if (
  !file.exists(file.path(
    output_directory, 
    paste0("snps_vep_v", ensembl_version, ".0_", genome_assembly, ".txt.gz")
  ))
) {
  input_docker <- gsub(
    pattern = ".*/outputs", 
    replacement = server_directory, 
    x = file.path(output_directory, "snps_locations.txt.gz")
  )
  output_docker <- paste0("snps_vep_", ensembl_version, ".0_", genome_assembly, ".txt")
  
  if (
    !file.exists(file.path(
      vep_cache[["docker"]], "homo_sapiens", paste0(ensembl_version, "_", genome_assembly)
    ))
  ) {
    system(paste0(
      "cd ", vep_cache[["docker"]],
      paste0(
        "curl -O ftp://ftp.ensembl.org/pub/release-", ensembl_version, 
        "/variation/vep/", ensembl_species, "_vep_", ensembl_version, "_", genome_assembly, ".tar.gz",
      ),
      paste0("tar xzf", ensembl_species, "_vep_", ensembl_version, "_", genome_assembly, ".tar.gz")
    ))
  }
  
  cat(paste(
    '#!/bin/bash',
    '\n\nchmod 777 -R',  gsub(".*/outputs", server_directory, output_directory),
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
    '--fork 100',
    '--force_overwrite',
    '--assembly', genome_assembly, 
    '--check_existing',
    '--no_check_alleles',
    '--symbol',
    '--output_file', file.path("/data_dir", output_docker),
    '&& cut -f 1-4,13-14', file.path("/data_dir", output_docker), 
    '| bgzip --thread 100 -f >', file.path("/data_dir", paste0(output_docker, ".gz")),
    '"',
    '\n\nchmod 775 -R',  gsub(".*/outputs", server_directory, output_directory),
    '\n'
  ), file = file.path(output_directory, "run_docker_vep.sh"))
} else {
  unlink(file.path(output_directory, paste0("snps_vep_v", ensembl_version, ".0_", genome_assembly, ".txt")))
}


### VEP ============================================================================================
if (
  !file.exists(file.path(
    output_directory, paste0("snps_vep_v", ensembl_version, ".0_", genome_assembly, "_formated.txt.gz")
  )) &
    file.exists(file.path(
      output_directory, paste0("snps_vep_v", ensembl_version, ".0_", genome_assembly, ".txt.gz")
    ))
) {
  vep_annotation <- fread(
    file = file.path(
      output_directory, paste0("snps_vep_v", ensembl_version, ".0_", genome_assembly, ".txt.gz")
    ), 
    skip = "#U"
  )[ 
    j = c("CHR", "POS") := tstrsplit(Location, ":", fixed = TRUE)
  ][
    j = (c("Gene", "Symbol", "rsid")) :=
      list(
        paste(unique(Gene), collapse = ";"),
        fifelse(
          test = grepl("SYMBOL=", Extra), 
          yes = paste(unique(gsub("^.*SYMBOL=([^;]*);.*$", "\\1", Extra)), collapse = ";"),
          no = NA_character_
        ),
        paste(unique(Existing_variation), collapse = ";")
      ), 
    by = "#Uploaded_variation"
  ][, .(CHR, POS, `#Uploaded_variation`, Gene, Symbol, rsid)]
  setnames(vep_annotation, "#Uploaded_variation", "chr_pos_ref_alt")
  fwrite(
    x = unique(vep_annotation), 
    file = file.path(
      output_directory, 
      paste0("snps_vep_v", ensembl_version, ".0_", genome_assembly, "_formated.txt.gz")
    )
  )
}


### Archive ========================================================================================
# if (!Sys.getenv("USER") %in% c("root", "") && file.exists("~/.fex/id")) {
#   local({
#     owd <- getwd()
#     setwd(normalizePath(output_directory))
#     archive_name <- file.path(
#       normalizePath(output_directory),
#       paste0(
#         format(Sys.Date(), format = "%Y%m%d"), "_",
#         project_name, "_", 
#         gsub("[0-9]+\\-", "", basename(output_directory)), ".zip"
#       )
#     )
#     zip(archive_name, files = list.files())
#     fex_out <- system(paste("fexsend", archive_name, "."), intern = TRUE)
#     unlink(archive_name)
#     setwd(owd)
#   })
# }


### Set chmod ======================================================================================
Sys.chmod(
  list.files(output_directory, full.names = TRUE),
  mode = "0775", use_umask = FALSE
)
Sys.chmod(
  list.files(output_directory, full.names = TRUE, recursive = TRUE, all.files = TRUE),
  mode = "0775", use_umask = FALSE
)
invisible(system(paste("chgrp -R staff", output_directory), intern = TRUE))


### Complete =======================================================================================
message("Success!", appendLF = TRUE)

