### global libraries ===============================================================================
library(targets)
library(tarchetypes)
library(here)
library(data.table)
# library(future)
# library(future.callr)

# targets::tar_renv(extras = "visNetwork", path = "scripts/_dependencies.R")


### project setup ==================================================================================
# Functions/scripts required: tar-qc_plink.R, tar-pval_trans.R
invisible(sapply(
  X = list.files(here("scripts", "tar-utils"), pattern = "^tar-.*R$", full.names = TRUE),
  FUN = source, echo = FALSE
))

# plan(future.callr::callr, workers = 40)
# plan(multicore, workers = 40)
# message(sprintf("Number of workers: %d", future::nbrOfWorkers()))
# setDTthreads(threads = 1)


### targets setup ==================================================================================
tar_setup <- list( # Setup project
  tar_target(project, sub("(.*)_[^_]*\\.Rproj$", "\\1", list.files(here(), pattern = ".Rproj$")), packages = "here"),
  tar_target(author, "Mickaël CANOUIL, *Ph.D.*"),
  tar_target(qc_data_path, "/disks/DATA/Projects"),
  tar_target(panel_data, "/disks/DATA/ExternalData"),
  tar_target(ga_run, "/disks/RUN/Array/Dmap/plink/PLINK_123456_0789/omni") # basename
)


### targets ========================================================================================
tar_genotype <- list( # Genotype Array (ga)
  tar_target(ga_imputation_panel,
    # command = file.path(panel_data, "1kg/hg19/Reference_genome/1000GP_Phase3_combined.legend"),
    command = file.path(panel_data, "HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab"),
    format = "file"
  ),
  tar_target(ga_ref1kg_panel,
    command = file.path(panel_data, "1kg/samples_description/integrated_call_samples_v3.20130502.ALL.panel"),
    format = "file"
  ),
  tar_target(ga_ref1kg_fasta,
    command = file.path(panel_data, "1kg/hg19/human_g1k_v37.fasta"),
    format = "file"
  ),
  tar_target(ga_ref1kg_population,
    command = file.path(panel_data, "1kg/samples_description/1kg_pop_description.tsv"),
    format = "file"
  ),
  tar_target(ga_params,
    command = list(
      input_files = ga_run,
      output_directory = here("outputs", "qc_ga"),
      array = "Illumina Omni2.5",
      callrate_samples = 0.95,
      callrate_snps = 0.95,
      sex_threshold = c(0.2, 0.8),
      heterozygosity_threshold = 4,
      maf_threshold = 0.01,
      hwe_pvalue = 0.0001,
      mendelian_samples = 0.05,
      mendelian_snps = 0.1,
      ibd_threshold = 0.2,
      population = NULL,
      pca_components = 10,
      pca_threshold = 3,
      max_labels = 15,
      check_bim_script = "https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip",
      imputation_ref = "HRC",
      imputation_panel = ga_imputation_panel,
      ref1kg_panel = ga_ref1kg_panel,
      ref1kg_fasta = ga_ref1kg_fasta,
      ref1kg_population = ga_ref1kg_population,
      ref1kg_genotypes = file.path(panel_data, "1kg/hg19/Genotypes/RawPLINK/ALL/")
    ),
    packages = "here"
  ),
  tar_target(ga_setup,
    command = create_ga_directory(ga_params[["output_directory"]])
  ),
  tar_target(ga_plink,
    command = download_plink(
      url = "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip",
      output_directory = ga_setup
    ),
    packages = "utils",
    format = "file"
  ),
  tar_target(ga_bcftools,
    command = "/usr/bin/bcftools",
    format = "file"
  ),
  tar_target(ga_vcftools,
    command = "/usr/bin/vcftools",
    format = "file"
  ),
  tar_target(ga_perl_imputation_check,
    command = download_perl_preimputation_check(
      url = ga_params[["check_bim_script"]],
      output_directory = ga_setup
    ),
    packages = "utils",
    format = "file"
  ),
  # tar_target(ga_plink_thread, # multi-threading
  #   command = sprintf("%s --threads %d", ga_plink, 20L)
  # ),
  tar_target(ga_bed,
    command = make_ga_bed(
      input = ga_params[["input_files"]],
      output = ga_setup,
      plink = ga_plink
    ),
    packages = "data.table",
    format = "file"
  ),
  tar_target(ga_bfile,
    command = sub("\\.bed$", "", grep("\\.bed$", ga_bed, value = TRUE))
  ),
  tar_target(ga_fam_data,
    command = read_fam(path = grep("\\.fam$", ga_bed, value = TRUE), project = project),
    packages = "data.table"
  ),
  tar_target(ga_samples_duplicated,
    command = count_duplicated_samples(ga_fam_data)
  ),
  tar_target(ga_callrate_samples,
    command = compute_callrate_ind(
      bfile = ga_bfile,
      callrate = ga_params[["callrate_samples"]],
      plink = ga_plink
    ),
    packages = "data.table"
  ),
  tar_target(ga_callrate_samples_plot,
    command = plot_callrate(
      data = ga_callrate_samples,
      callrate = ga_params[["callrate_samples"]],
      max_labels = ga_params[["max_labels"]],
      type = "sample"
    ),
    packages = c("data.table", "ggplot2", "ggrepel", "scales")
  ),
  tar_target(ga_callrate_samples_leq_geq,
    command = compute_snps_maf_leq_geq(
      bfile = ga_bfile,
      maf = ga_params[["maf_threshold"]],
      callrate = ga_params[["callrate_samples"]],
      plink = ga_plink
    ),
    packages = "data.table"
  ),
  tar_target(ga_check_sex,
    command = check_genotypic_sex(
      bfile = ga_bfile,
      callrate_data = ga_callrate_samples,
      fam_data = ga_fam_data,
      sex_threshold = ga_params[["sex_threshold"]],
      plink = ga_plink
    ),
    packages = "data.table"
  ),
  tar_target(ga_check_sex_plot,
    command = plot_check_genotypic_sex(
      data = ga_check_sex,
      callrate = ga_params[["callrate_samples"]],
      sex_threshold = ga_params[["sex_threshold"]],
      max_labels = ga_params[["max_labels"]]
    ),
    packages = c("data.table", "ggplot2", "ggrepel", "scales")
  ),
  tar_target(ga_het_samples_all,
    command = compute_snps_het(
      bfile = ga_bfile,
      heterozygosity_threshold = ga_params[["heterozygosity_threshold"]],
      callrate = ga_params[["callrate_samples"]],
      callrate_data = ga_callrate_samples,
      plink = ga_plink,
      snps = NULL,
      what = NULL
    ),
    packages = "data.table"
  ),
  tar_target(ga_het_samples_leq_geq,
    command = compute_snps_het_leq_geq(
      bfile = ga_bfile,
      heterozygosity_threshold = ga_params[["heterozygosity_threshold"]],
      maf = ga_params[["maf_threshold"]],
      callrate = ga_params[["callrate_samples"]],
      callrate_data = ga_callrate_samples,
      plink = ga_plink
    ),
    packages = "data.table"
  ),
  tar_target(ga_het_samples,
    command = compute_snps_het_all(ga_het_samples_all, ga_het_samples_leq_geq),
    packages = "data.table"
  ),
  tar_target(ga_het_samples_plot,
    command = plot_het_ind(
      data = ga_het_samples,
      heterozygosity_threshold = ga_params[["heterozygosity_threshold"]],
      callrate = ga_params[["callrate_samples"]],
      maf = ga_params[["maf_threshold"]],
      max_labels = ga_params[["max_labels"]]
    ),
    packages = c("data.table", "ggplot2", "ggrepel", "scales")
  ),
  tar_target(ga_relatedness,
    command = compute_relatedness(
      bfile = ga_bfile,
      maf = ga_params[["maf_threshold"]],
      ibd = ga_params[["ibd_threshold"]],
      plink = ga_plink
    ),
    packages = "data.table"
  ),
  tar_target(ga_related_samples,
    command = compute_related_samples_tab(
      relatedness = ga_relatedness,
      callrate_samples = ga_callrate_samples
    ),
    packages = "data.table"
  ),
  tar_target(ga_pca_ethnicity,
    command = compute_pca_ethnicity(
      bfile = ga_bfile,
      maf = ga_params[["maf_threshold"]],
      hwe = ga_params[["hwe_pvalue"]],
      ref1kg_genotypes = ga_params[["ref1kg_genotypes"]],
      plink = ga_plink
    ),
    packages = c("data.table", "future", "future.apply", "flashpcaR")
  ),
  tar_target(ga_pca_ethnicity_tidy,
    command = tidy_pca_ethnicity(
      data = ga_pca_ethnicity[["vectors"]],
      ref1kg_panel = ga_params[["ref1kg_panel"]]
    ) ,
    packages = "data.table"
  ),
  tar_target(ga_pca_ethnicity_tidy_plot,
    command = plot_pca_ethnicty(
      data = ga_pca_ethnicity_tidy,
      pve = ga_pca_ethnicity[["pve"]],
      loadings = ga_pca_ethnicity[["loadings"]]
    ),
    packages = c("data.table", "ggplot2", "ggforce", "concaveman", "ggtext")
  ),
  tar_target(ga_samples_exclude,
    command = compute_samples_to_exclude(
      callrate_data = ga_callrate_samples,
      sexcheck_data = ga_check_sex,
      heterozygosity_data = ga_het_samples,
      relatedness_data = ga_relatedness
    ),
    packages = "data.table"
  ),
  tar_target(ga_bed_good_samples,
    command = compute_bed_good_samples(
      bfile = ga_bfile,
      exclude = ga_samples_exclude[Status %in% "Exclude"],
      temp_directory = ga_setup,
      plink = ga_plink
    ),
    packages = "data.table",
    format = "file"
  ),
  tar_target(ga_bfile_good_samples,
    command = sub("\\.bed$", "", grep("\\.bed$", ga_bed_good_samples, value = TRUE))
  ),
  tar_target(ga_callrate_snp,
    command = compute_callrate_snp(
      bfile = ga_bfile_good_samples,
      callrate = ga_params[["callrate_snps"]],
      plink = ga_plink
    ),
    packages = "data.table"
  ),
  tar_target(ga_callrate_snp_plot,
    command = plot_callrate(
      data = ga_callrate_snp,
      callrate = ga_params[["callrate_snps"]],
      max_labels = ga_params[["max_labels"]],
      type = "snp"
    ),
    packages = c("data.table", "ggplot2", "ggrepel", "scales")
  ),
  tar_target(ga_duplicated_snps,
    command = compute_duplicated_snp(
      bfile = ga_bfile_good_samples
    ),
    packages = "data.table"
  ),
  tar_target(ga_hwe_snp,
    command = compute_hwe_snp(
      bfile = ga_bfile_good_samples,
      hwe = ga_params[["hwe_pvalue"]],
      plink = ga_plink
    ),
    packages = "data.table"
  ),
  tar_target(ga_hwe_snp_plot,
    command = plot_hwe_snp(
      data = ga_hwe_snp,
      hwe = ga_params[["hwe_pvalue"]]
    ),
    packages = c("data.table", "ggplot2", "scales")
  ),
  tar_target(ga_maf_snp,
    command = compute_maf_snp(
      bfile = ga_bfile_good_samples,
      maf = ga_params[["maf_threshold"]],
      plink = ga_plink
    ),
    packages = "data.table"
  ),
  tar_target(ga_maf_snp_plot,
    command = plot_maf_snp(data = ga_maf_snp),
    packages = c("data.table", "ggplot2", "scales")
  ),
  tar_target(ga_bed_good_samples_variants,
    command = compute_bed_good_variants(
      bfile = ga_bfile_good_samples,
      exclude = list(
        ga_callrate_snp,
        ga_duplicated_snps,
        ga_hwe_snp
      ),
      temp_directory = ga_setup,
      project = project,
      plink = ga_plink
    ),
    packages = "data.table",
    format = "file"
  ),
  tar_target(ga_bfile_good_samples_variants,
    command = sub("\\.bed$", "", grep("\\.bed$", ga_bed_good_samples_variants, value = TRUE))
  ),
  tar_target(ga_vcf_imputation,
    command = compute_vcf_imputation(
      bfile = ga_bfile_good_samples_variants,
      ref = ga_params[["imputation_ref"]],
      ref_panel = ga_params[["imputation_panel"]],
      ref1kg_fasta = ga_params[["ref1kg_fasta"]],
      perl_script = ga_perl_imputation_check,
      temp_directory = ga_setup,
      project = project,
      plink = ga_plink,
      bcftools = ga_bcftools
    ),
    packages = "data.table",
    format = "file"
  ),
  tar_target(ga_vcf_imputation_dim,
    command = compute_vcf_dim(grep("\\.vcf.gz$", ga_vcf_imputation, value = TRUE)),
    packages = c("data.table", "R.utils")
  ),
  tar_target(ga_export_directory,
    command = create_ga_export_directory(
      path = qc_data_path,
      project = project,
      array = sub(".* ", "", ga_params[["array"]])
    )
  ),
  tar_target(ga_imputed_vcf,
    command = list_imputed_vcf(ga_export_directory),
    format = "file",
    error = "continue"
  ),
  tar_target(ga_vcf_imputed_uptodate,
    command = is_vcf_imputed_uptodate(
      pre = ga_vcf_imputation,
      post = ga_imputed_vcf
    )
  ),
  tar_target(ga_vcf_imputed_qc,
    command = compute_vcf_imputed_qc(
      vcf = ga_imputed_vcf,
      vcftools = ga_vcftools,
      uptodate = ga_vcf_imputed_uptodate
    ),
    packages = c("data.table", "future.apply"),
    error = "continue"
  ),
  tar_target(ga_vcf_imputed_qc_info,
    command = compute_vcf_imputed_qc_info(
      data = ga_vcf_imputed_qc,
      uptodate = ga_vcf_imputed_uptodate
    ),
    packages = c("data.table", "scales"),
    error = "continue"
  ),
  tar_target(ga_vcf_imputed_qc_af,
    command = compute_vcf_imputed_qc_af(
      data = ga_vcf_imputed_qc,
      uptodate = ga_vcf_imputed_uptodate
    ),
    packages = c("data.table", "scales"),
    error = "continue"
  ),
  tar_target(ga_vcf_imputed_qc_plot_export,
    command = save_plot_vcf_imputed_qc(
      data = ga_vcf_imputed_qc,
      path = ga_setup,
      uptodate = ga_vcf_imputed_uptodate
    ),
    packages = c(
      "data.table", "ggplot2", "patchwork", "scales",
      "utils", "ragg"
    ),
    format = "file",
    error = "continue"
  ),
  tar_render(ga_qc_report,
    path = here("scripts/slides/genotype_array_qc_report.Rmd"),
    output_dir = here("reports"),
    params = ga_params,
    packages = c(
      "xaringan",
      "here", "knitr", "ragg", "ggplot2", "ggtext",
      "patchwork", "data.table", "gt", "scales",
      "targets"
    )
  ),
  tar_target(ga_save_qc,
    command = save_ga_qc_data(
      from = dirname(ga_setup),
      report = here(grep("\\.html$", ga_qc_report, value = TRUE)),
      imputation = ga_vcf_imputed_qc_plot_export,
      exclusion_check = ga_samples_exclude,
      relatedness = ga_related_samples,
      ethnicity = ga_pca_ethnicity_tidy,
      to = ga_export_directory
    )
  )
)

list(
  tar_setup,
  tar_genotype
)
