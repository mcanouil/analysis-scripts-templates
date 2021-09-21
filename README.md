
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Mickaël Canouil’s Analysis Scripts Templates

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- [![GitHub
tag](https://img.shields.io/github/tag/mcanouil/analysis_templates.svg?label=latest%20tag&include_prereleases)](https://github.com/mcanouil/analysis_templates) -->
<!-- badges: end -->

Script templates to work with omics datasets.

## R

  - [default.R](R/default.R)

### R/functions

  - [tar-crossmap.R](R/functions/tar-crossmap.R)
  - [tar-ewas.R](R/functions/tar-ewas.R)
  - [tar-gwas.R](R/functions/tar-gwas.R)
  - [tar-manhattan.R](R/functions/tar-manhattan.R)
  - [tar-plink2.R](R/functions/tar-plink2.R)
  - [tar-pval\_trans.R](R/functions/tar-pval_trans.R)
  - [tar-qc\_idats.R](R/functions/tar-qc_idats.R)
  - [tar-qc\_plink.R](R/functions/tar-qc_plink.R)
  - [tar-sub\_chunk.R](R/functions/tar-sub_chunk.R)
  - [tar-vep.R](R/functions/tar-vep.R)

### R/scripts

  - [01-design.R](R/scripts/01-design.R)
  - [02-qc\_idats\_snps.Rmd](R/scripts/02-qc_idats_snps.Rmd)
  - [02-qc\_idats.Rmd](R/scripts/02-qc_idats.Rmd)
  - [03-qc\_plink.Rmd](R/scripts/03-qc_plink.Rmd)
  - [04-qc\_impute.Rmd](R/scripts/04-qc_impute.Rmd)
  - [05-ethnicity.R](R/scripts/05-ethnicity.R)
  - [06-omni\_to\_hg38.R](R/scripts/06-omni_to_hg38.R)
  - [07-save\_outputs.R](R/scripts/07-save_outputs.R)
  - [08-vep\_vcf\_docker.R](R/scripts/08-vep_vcf_docker.R)
  - [09-gwas.R](R/scripts/09-gwas.R)
  - [10-ewas.R](R/scripts/10-ewas.R)
  - [11-twas.R](R/scripts/11-twas.R)
  - [12-ora.R](R/scripts/12-ora.R)
  - [13-gsea.R](R/scripts/13-gsea.R)
  - [14-mqtl.R](R/scripts/14-mqtl.R)
  - [15-eqtl.R](R/scripts/15-eqtl.R)
  - [16-eqtm\_estimated\_time.R](R/scripts/16-eqtm_estimated_time.R)
  - [16-eqtm.R](R/scripts/16-eqtm.R)

### R/targets

  - [tar\_crossmap.R](R/targets/tar_crossmap.R)
  - [tar\_ewas.R](R/targets/tar_ewas.R)
  - [tar\_genotype.R](R/targets/tar_genotype.R)
  - [tar\_gsea.R](R/targets/tar_gsea.R)
  - [tar\_gwas.R](R/targets/tar_gwas.R)
  - [tar\_methylation.R](R/targets/tar_methylation.R)
  - [tar\_ora.R](R/targets/tar_ora.R)
  - [tar\_vep.R](R/targets/tar_vep.R)

## Rmarkdown

Note: `assets` directory from
[github.com/umr1283/xaringan-template](https://github.com/umr1283/xaringan-template)
is required\!

  - [ewas\_report.Rmd](Rmarkdown/ewas_report.Rmd)
  - [genotype\_array\_qc\_report.Rmd](Rmarkdown/genotype_array_qc_report.Rmd)
  - [gwas\_report.Rmd](Rmarkdown/gwas_report.Rmd)
  - [methylation\_array\_qc\_report.Rmd](Rmarkdown/methylation_array_qc_report.Rmd)
