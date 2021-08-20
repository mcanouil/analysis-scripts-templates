message(timestamp(quiet = TRUE))
### Project Setup ==================================================================================
library(here)
project_name <- sub("(.*)_*\\..*", "\\1", list.files(here(), pattern = ".Rproj$"))
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


### Functions ======================================================================================
source("functions/tar-vep.R") 


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
output <- get_symbol_vep(
  input = file.path(output_directory, "snps_locations.txt.gz"),
  output_directory = here::here("outputs"),
  genome_assembly = genome_assembly,
  ensembl_version = ensembl_version,
  ensembl_species = ensembl_species,
  vep_cache = vep_cache
)
message(sprintf(
  "Run the following command before: 'nohup bash %s/run_docker_vep.sh > /dev/null &'", 
  output_directory
))

format_symbol_vep(output)


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
message(timestamp(quiet = TRUE))

