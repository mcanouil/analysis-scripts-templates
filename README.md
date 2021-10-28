
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Mickaël Canouil’s Analysis Scripts Templates

<!-- badges: start -->

[![License](https://img.shields.io/github/license/mcanouil/analysis-scripts-templates)](LICENSE)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)]()
<!-- [![GitHub
tag](https://img.shields.io/github/tag/mcanouil/analysis-scripts-templates.svg?label=latest%20tag&include_prereleases)](https://github.com/mcanouil/analysis-scripts-templates) -->
<!-- badges: end -->

Script templates to work with omics datasets.

## R

-   [default.R](R/default.R)

### R/functions

-   [tar-crossmap.R](R/functions/tar-crossmap.R)
-   [tar-ewas.R](R/functions/tar-ewas.R)
-   [tar-gwas.R](R/functions/tar-gwas.R)
-   [tar-manhattan.R](R/functions/tar-manhattan.R)
-   [tar-meqtl.R](R/functions/tar-meqtl.R)
-   [tar-plink2.R](R/functions/tar-plink2.R)
-   [tar-pval_trans.R](R/functions/tar-pval_trans.R)
-   [tar-qc_idats.R](R/functions/tar-qc_idats.R)
-   [tar-qc_plink.R](R/functions/tar-qc_plink.R)
-   [tar-sub_chunk.R](R/functions/tar-sub_chunk.R)
-   [tar-vep.R](R/functions/tar-vep.R)

### R/scripts

-   [01-design.R](R/scripts/01-design.R)
-   [02-qc_idats.Rmd](R/scripts/02-qc_idats.Rmd)
-   [02-qc_idats_snps.Rmd](R/scripts/02-qc_idats_snps.Rmd)
-   [03-qc_plink.Rmd](R/scripts/03-qc_plink.Rmd)
-   [04-qc_impute.Rmd](R/scripts/04-qc_impute.Rmd)
-   [05-ethnicity.R](R/scripts/05-ethnicity.R)
-   [06-omni_to_hg38.R](R/scripts/06-omni_to_hg38.R)
-   [07-save_outputs.R](R/scripts/07-save_outputs.R)
-   [08-vep_vcf_docker.R](R/scripts/08-vep_vcf_docker.R)
-   [09-gwas.R](R/scripts/09-gwas.R)
-   [10-ewas.R](R/scripts/10-ewas.R)
-   [11-twas.R](R/scripts/11-twas.R)
-   [12-ora.R](R/scripts/12-ora.R)
-   [13-gsea.R](R/scripts/13-gsea.R)
-   [14-meqtl.R](R/scripts/14-meqtl.R)
-   [15-eqtl.R](R/scripts/15-eqtl.R)
-   [16-eqtm.R](R/scripts/16-eqtm.R)
-   [16-eqtm_estimated_time.R](R/scripts/16-eqtm_estimated_time.R)

### R/targets

-   [tar_crossmap.R](R/targets/tar_crossmap.R)
-   [tar_ewas.R](R/targets/tar_ewas.R)
-   [tar_genotype.R](R/targets/tar_genotype.R)
-   [tar_gsea.R](R/targets/tar_gsea.R)
-   [tar_gwas.R](R/targets/tar_gwas.R)
-   [tar_meqtl.R](R/targets/tar_meqtl.R)
-   [tar_methylation.R](R/targets/tar_methylation.R)
-   [tar_ora.R](R/targets/tar_ora.R)
-   [tar_vep.R](R/targets/tar_vep.R)

## Rmarkdown

Note: `assets` directory from
[github.com/umr1283/xaringan-template](https://github.com/umr1283/xaringan-template)
is required!

-   [ewas_report.Rmd](Rmarkdown/ewas_report.Rmd)
-   [genotype_array_qc_report.Rmd](Rmarkdown/genotype_array_qc_report.Rmd)
-   [gwas_report.Rmd](Rmarkdown/gwas_report.Rmd)
-   [methylation_array_qc_report.Rmd](Rmarkdown/methylation_array_qc_report.Rmd)
