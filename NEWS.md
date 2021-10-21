# analysis_templates (development version)

* `R/targets/`
  + `tar_crossmap.R`, `targets` setup to upgrade or downgrade genome assembly of VCF files.
  + `tar_ewas.R`, `targets` setup to perform multiple EWAS.
  + `tar_genotype.R`, `targets` setup to compute PLINK files quality-control.
  + `tar_gsea.R`, `targets` setup to compute gene-set enrichment analysis.
  + `tar_gwas.R`, `targets` setup to perform multiple GWAS using PLINK2 software.
  + `tar_methylation.R`, `targets` setup to compute idats (methylation) files quality-control.
  + `tar_ora.R`, `targets` setup to compute over-representation analysis.
  + `tar_vep.R`, `targets` setup to compute VEP annotations.

* `R/functions/`
  + `tar-crossmap.R`, functions, used in `tar_crossmap.R`, to upgrade or downgrade genome assembly of VCF files.
  + `tar-ewas.R`, function, used in `tar_ewas.R`, to perform multiple EWAS.
  + `tar-gwas.R`, function, used in `tar_gwas.R`, to perform multiple GWAS using PLINK2 software.
  + `tar-manhattan.R`, functions to draw Manhattan plot using `ggplot2`.
  + `tar-plink2.R`n function to download PLINK2 binary.
  + `tar-pval_trans.R`, p-value `ggplot2` axis transformation.
  + `tar-qc_idats.R`, functions, used in `tar_methylation.R`, to compute idats (methylation) files quality-control.
  + `tar-qc_plink.R`, functions, used in `tar_genotype.R`, to compute PLINK files quality-control.
  + `tar-sub_chunk.R`, function to create sub-chunk inside a chunk in RMarkdown.
  + `tar-vep.R`, functions, used in `tar_vep.R`, to prepare Docker command and to format VEP output.

* `R/rmarkdown/` (`assets` directory from [github.com/umr1283/xaringan-template](https://github.com/umr1283/xaringan-template) is required)
  + `ewas_report.Rmd`, `xaringan` report template associated with `R/targets/tar_ewas.R` and `R/functions/tar-ewas.R`.
  + `genotype_array_qc_report.Rmd`, `xaringan` report template associated with `R/targets/tar_genotype.R` and `R/functions/tar-qc_plink.R`.
  + `gwas_report.Rmd`, `xaringan` report template associated with `R/targets/tar_gwas.R` and `R/functions/tar-gwas.R`.
  + `methylation_array_qc_report.Rmd`, `xaringan` report template associated with `R/targets/tar_methylation.R` and `R/functions/tar-qc_idats.R`.

* `R/scripts`
  + `**-****.R`,
    - fix project name to not include file extension.
    - trim trailing spaces.
  + `default.R`,
    - default R script.
  + `01-design.R`,
    - skeleton script make batch design for omics (currently "empty").
  + `02-qc_idats.Rmd`,
    - YAML header for Rmarkdown template in `umr1283/dmapaq`.
  + `02-qc_idats_snps.Rmd`,
    - YAML header for Rmarkdown template in `umr1283/dmapaq` (mostly for mQTL analysis).
  + `03-qc_plink.Rmd`,
    - skeleton script or YAML header for Rmarkdown template in `umr1283/dgapaq`.
  + `04-qc_impute.Rmd`,
    - skeleton script or YAML header for Rmarkdown template in `umr1283/dgapaq`.
  + `05-ethnicity.R`,
    - skeleton script to perform ethnicity inference using 1,000 Genomes and HRC imputed VCF files.
  + `06-omni_to_hg38.R`,
    - skeleton R script to upgrade genotyping arrays to GRCh38 genome assembly.
  + `07-save_outputs.R`,
    - skeleton R script to store "sensitive"/"fragile" data to a more secure disk.
  + `08-vep_vcf_docker.R`,
    - skeleton R script with iterative steps to retrieve RSID and gene's symbol using VEP (Docker).
    - add function to prepare docker command and to format output.
  + `09-gwas.R`,
    - skeleton R script to perform a GWAS using `PLINK2`  with multiple traits.
  + `10-ewas.R`,
    - skeleton R script to perform a EWAS using `limma` with multiple traits.
    - tweak & refactor some code.
    - use as.matrix from `data.table`.
  + `11-twas.R`,
    - skeleton R script to perform a GWAS using `DESeq2` with multiple traits.
    - fix #2 (@Ning-L).
    - ensure Ensembl ID does not contain version.
  + `12-ora.R`,
    - add gene symbols for KEGG.
  + `13-gsea.R` and `R/tar_gsea.R`,
    - add genes set and peripheral sets (genes set - core set).
    - add gene symbols in addition to Ensembl, Uniprot or Entrez gene IDs.
    - fix ties using pvalue as secondary parameter in sorting and removing `.Machine$double.eps` for the second value.
  + `14-mqtl.R`,
    - skeleton R script to perform mQTL analysis (methylation Vs. genotypes).
    - fix missing variables.
    - tweak & refactor some code.
    - use as.matrix from `data.table`.
  + `15-eqtl.R`,
    - keleton R script to perform eQTL analysis (expression Vs. genotypes).
    - fix missing variables.
  + `16-eqtm.R`,
    - new expression-methylation analysis.
    - tweak & refactor some code.
    - isoform/transcript or gene level.
  + `16-eqtm_estimated_time.R`, new script to estimate computation time from `16-eqtm.R`'s log.
