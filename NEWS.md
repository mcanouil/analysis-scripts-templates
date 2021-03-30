# analysis_templates (development version)

+ `R/13-mqtl.R`, 
    + Skeleton R script to perform mQTL analysis (methylation Vs. genotypes).
+ `R/14-eqtl.R`, 
    + Skeleton R script to perform eQTL analysis (expression Vs. genotypes.
+ `R/12-ora_gsea.R`, 
    + Add genes set and peripheral sets (genes set - core set).
    + Fix typos.

# analysis_templates 0.1.0

+ `R/default.R`, default R script.
+ `R/01-design.R`, skeleton script make batch design for omics (currently "empty").
+ `R/02-qc_idats.Rmd`, YAML header for Rmarkdown template in `umr1283/dmapaq`.
+ `R/02-qc_idats_snps.Rmd`, YAML header for Rmarkdown template in `umr1283/dmapaq` (mostly for mQTL analysis).
+ `R/03-qc_plink.Rmd`, skeleton script or YAML header for Rmarkdown template in `umr1283/dgapaq`.
+ `R/04-qc_impute.Rmd`, skeleton script or YAML header for Rmarkdown template in `umr1283/dgapaq`.
+ `R/05-ethnicity.R`, skeleton script to perform etnicity inference using 1,000 Genomes and HRC imputed VCF files.
+ `R/06-omni_to_hg38.R`, skeleton R script to upgrade genotyping arrays to GRCh38 genome assembly.
+ `R/07-save_outputs.R`, skeleton R script to store "sensitive"/"fragile" data to a more secure disk.
+ `R/08-vep_vcf_docker.R`, skeleton R script with iterative steps to retrieve RSID and gene's symbol using VEP (Docker).
+ `R/09-gwas.R`, skeleton R script to perform a GWAS using `PLINK2`  with multiple traits.
+ `R/10-ewas.R`, skeleton R script to perform a EWAS using `limma` with multiple traits.
+ `R/11-twas.R`, skeleton R script to perform a GWAS using `DESeq2` with multiple traits.
+ `R/12-ora_gsea.R`, skeleton R script to perform enrichment and gene-set enrichment analyses using `clusterProfiler `.
