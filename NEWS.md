# analysis_templates (development version)

+ `target/`
    + `tar_ora.R`, `targets` setup to compute over-representation analysis.
    + `tar_gsea.R`, `targets` setup to compute gene-set enrichment analysis.
    + `tar_methylation.R`, `targets` setup to compute idats (methylation) files quality-control.
    + `tar-qc_idats`, functions, used in `tar_methylation.R`, to compute idats (methylation) files quality-control.
    + `tar_genotype.R`, `targets` setup to compute PLINK files quality-control.
    + `tar-qc_plink`, functions, used in `tar_genotype.R`, to compute PLINK files quality-control.
    + `tar-pval_trans`, p-value ggplot2 axis transformation.
    + `tar-sub_chunk`, function to create sub-chunk inside a chunk in RMarkdown.

+ `R/`
    + `**-****.R`, fix project name to not include file extension. 
    + `10-ewas.R`, tweak & refactor some code.
    + `11-twas.R`, 
        - fix #2 (@Ning-L).
        - ensure Ensembl ID does not contain version.
    + `14-mqtl.R`, 
        - fix missing variables.
        - tweak & refactor some code.
    + `15-eqtl.R`, fix missing variables.
    + `16-eqtm.R`, 
        - new expression-methylation analysis.
        - tweak & refactor some code.
    + `16-eqtm_estimated_time.R`, new script to estimate computation time from `16-eqtm.R`'s log.

# analysis_templates 0.2.0

+ `R/14-mqtl.R`, 
    + Skeleton R script to perform mQTL analysis (methylation Vs. genotypes).
+ `R/15-eqtl.R`, 
    + Skeleton R script to perform eQTL analysis (expression Vs. genotypes).
+ `R/12-ora.R` and `R/tar_ora.R`, 
    + Add gene symbols for KEGG.
    + Add ORA as a simple R script or `targets` setup.
+ `R/13-gsea.R` and `R/tar_gsea.R`, 
    + Add genes set and peripheral sets (genes set - core set).
    + Add gene symbols in addition to Ensembl, Uniprot or Entrez gene IDs.
    + Fix ties using pvalue as secondary parameter in sorting and removing `.Machine$double.eps` for the second value.
    + Add GSEA as a simple R script or `targets` setup.

# analysis_templates 0.1.0

+ `R/default.R`, default R script.
+ `R/01-design.R`, skeleton script make batch design for omics (currently "empty").
+ `R/02-qc_idats.Rmd`, YAML header for Rmarkdown template in `umr1283/dmapaq`.
+ `R/02-qc_idats_snps.Rmd`, YAML header for Rmarkdown template in `umr1283/dmapaq` (mostly for mQTL analysis).
+ `R/03-qc_plink.Rmd`, skeleton script or YAML header for Rmarkdown template in `umr1283/dgapaq`.
+ `R/04-qc_impute.Rmd`, skeleton script or YAML header for Rmarkdown template in `umr1283/dgapaq`.
+ `R/05-ethnicity.R`, skeleton script to perform ethnicity inference using 1,000 Genomes and HRC imputed VCF files.
+ `R/06-omni_to_hg38.R`, skeleton R script to upgrade genotyping arrays to GRCh38 genome assembly.
+ `R/07-save_outputs.R`, skeleton R script to store "sensitive"/"fragile" data to a more secure disk.
+ `R/08-vep_vcf_docker.R`, skeleton R script with iterative steps to retrieve RSID and gene's symbol using VEP (Docker).
+ `R/09-gwas.R`, skeleton R script to perform a GWAS using `PLINK2`  with multiple traits.
+ `R/10-ewas.R`, skeleton R script to perform a EWAS using `limma` with multiple traits.
+ `R/11-twas.R`, skeleton R script to perform a GWAS using `DESeq2` with multiple traits.
+ `R/12-ora_gsea.R`, skeleton R script to perform enrichment and gene-set enrichment analyses using `clusterProfiler `.
