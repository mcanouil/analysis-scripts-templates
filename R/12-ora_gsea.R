message(timestamp(quiet = TRUE))
### Project Setup ==================================================================================
library(here)
project_name <- gsub("(.*)_.*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "12-ora_gsea")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

organism <- c(
  "reactome" = "human",
  "go" = "org.Hs.eg.db",
  "kegg" = "hsa"
)

pvalue_gene_ora <- 0.05
fdr_pathway_ora <- 0.10

pvalue_gene_gsea <- 0.05
fdr_pathway_gsea <- 0.10


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(ragg)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
  library(data.table)
  library(clusterProfiler)
  library(ReactomePA)
  library(organism[["go"]], character.only = TRUE)
  library(glue)
  library(writexl)
})


### Tables and Figures Theme =======================================================================
options(
  ggplot2.discrete.colour = function(...) scale_colour_viridis_d(..., begin = 0.15, end = 0.85),
  ggplot2.discrete.fill = function(...) scale_fill_viridis_d(..., begin = 0.15, end = 0.85),
  ggplot2.continuous.colour = function(...) scale_colour_viridis_c(..., begin = 0.15, end = 0.85),
  ggplot2.continuous.fill = function(...) scale_fill_viridis_c(..., begin = 0.15, end = 0.85)
)
theme_set(
  theme_minimal(base_family = "Times") +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown(face = "italic", size = rel(0.8)),
      plot.caption = element_markdown(face = "italic", size = rel(0.5)),
      axis.title.x = element_markdown(),
      axis.text.x = element_markdown(),
      axis.title.y = element_markdown(),
      axis.text.y = element_markdown()
    )
)


### Analysis =======================================================================================
results_file <- here("outputs", "11-twas", paste0(project_name, "_DE_status_gene.csv.gz"))

## Over-representation -----------------------------------------------------------------------------
invisible(lapply(
  X = results_file,
  FUN = function(.file) {
    results <- fread(.file)[order(pvalue)]
    enrich_sets <- list(
      "Reactome" = enrichPathway(
        gene = unique(na.exclude(unlist(strsplit(results[pvalue < pvalue_gene_ora][["entrezgene_id"]], ";")))), 
        universe = unique(na.exclude(unlist(strsplit(results[["entrezgene_id"]], ";")))),
        organism = organism[["reactome"]],
        pvalueCutoff = fdr_pathway_ora, 
        pAdjustMethod = "BH",
        readable = TRUE
      ),
      "Gene Ontology Biological Process" = enrichGO(
        gene = unique(results[pvalue < pvalue_gene_ora][["ensembl_gene_id"]]),
        universe = unique(results[["ensembl_gene_id"]]),
        OrgDb = get(organism[["go"]]),
        keyType = "ENSEMBL",
        ont = "BP",
        pvalueCutoff = fdr_pathway_ora,
        pAdjustMethod = "BH",
        readable = TRUE
      ),
      "Gene Ontology Cellular Component" = enrichGO(
        gene = unique(results[pvalue < pvalue_gene_ora][["ensembl_gene_id"]]),
        universe = unique(results[["ensembl_gene_id"]]),
        OrgDb = get(organism[["go"]]),
        keyType = "ENSEMBL",
        ont = "CC",
        pvalueCutoff = fdr_pathway_ora,
        pAdjustMethod = "BH",
        readable = TRUE
      ),
      "Gene Ontology Molecular Function" = enrichGO(
        gene = unique(results[pvalue < pvalue_gene_ora][["ensembl_gene_id"]]),
        universe = unique(results[["ensembl_gene_id"]]),
        OrgDb = get(organism[["go"]]),
        keyType = "ENSEMBL",
        ont = "MF",
        pvalueCutoff = fdr_pathway_ora,
        pAdjustMethod = "BH",
        readable = TRUE
      ),
      "KEGG" = enrichKEGG(
        gene = unique(na.exclude(unlist(strsplit(results[pvalue < pvalue_gene_ora][["uniprotswissprot"]], ";")))), 
        universe = unique(na.exclude(unlist(strsplit(results[["uniprotswissprot"]], ";")))),
        organism = organism[["kegg"]], 
        keyType = "uniprot", 
        pvalueCutoff = fdr_pathway_ora, 
        pAdjustMethod = "BH"
      )
    )
    
    write_xlsx(
      x = setNames(lapply(enrich_sets, FUN = function(.enrich) {
        if (is.null(.enrich)) return(data.frame())
        .enrich@result
      }), gsub("Gene Ontology", "GO", names(enrich_sets))), 
      path = file.path(output_directory, "over_representation.xlsx")
    )
    
    enrich_sets <- enrich_sets[!sapply(enrich_sets, is.null)]
    
    if (
      all(
        sapply(
          X = enrich_sets, 
          FUN = function(.enrich) nrow(setDT(.enrich@result)[p.adjust < fdr_pathway_ora])
        ) == 0
      )
    ) {
      return(NULL)
    }
    
    plot_scales <- rbindlist(lapply(enrich_sets, FUN = function(.enrich) {
      setDT(.enrich@result)[
        p.adjust < fdr_pathway_ora
      ][
        j = list(
          p.adjust = max(round(-log10(c(fdr_pathway_ora, p.adjust)), digits = 0)),
          Count = max(c(Count, 0))
        )
      ]
    }))[j = lapply(.SD, max)] + 1
    
    enrich_plots <- lapply(enrich_sets, FUN = function(.enrich) {
      .dt <- setDT(.enrich@result)[
        p.adjust < fdr_pathway_ora
      ][
        j = GeneRatio := sapply(GeneRatio, function(.x) eval(parse(text = .x)))
      ][
        order(GeneRatio)
      ][
        j = Description := factor(Description, levels = unique(Description))
      ]
      
      if (nrow(.dt) == 0) {
        return(
          ggplot() +
            annotate(geom = "text", x = 1, y = 1, label = paste0("No enriched terms found\n(FDR < ", fdr_pathway_ora, ").")) +
            theme(
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              panel.grid = element_blank()
            )
        )
      }
      
      is_too_big <- nrow(.dt) > 25
      
      if (is_too_big) .dt <- first(.dt[order(pvalue)], 25)
      
      ggplot(data = .dt) +
        aes(x = GeneRatio, y = Description, colour = -log10(p.adjust), size = Count) +
        geom_point() +
        scale_y_discrete(
          labels = function(.x) {
            .x <- ifelse(
              nchar(.x) > 50,
              gsub("(.{0,50}).*", "\\1 [...]", .x),
              .x
            )
            gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", .x, perl = TRUE)
          }
        ) +
        scale_colour_continuous(limits = c(-log10(fdr_pathway_ora), plot_scales[["p.adjust"]])) +
        scale_size_continuous(limits = c(1, plot_scales[["Count"]]), range = c(0.1, 3)) +
        labs(
          x = if (is_too_big) "Gene Ratio (Top 25)" else "Gene Ratio", 
          y = NULL, 
          colour = "-log<sub>10</sub>(FDR)<br>(Benjamini-Hochberg)",
          size = "Number of Genes"
        ) +
        theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
        theme(
          plot.title.position = "plot",
          plot.caption.position = "plot",
          plot.title = element_markdown(),
          plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
          plot.caption = element_markdown(face = "italic", size = rel(0.75)),
          axis.title.x = element_markdown(),
          axis.text.x = element_markdown(),
          axis.title.y = element_markdown(),
          axis.text.y = element_markdown(size = rel(0.75)),
          panel.grid.minor = element_blank(),
          egend.title = element_markdown()
        ) +
        guides(
          colour = guide_colourbar(
            title.position = "top", 
            title.hjust = 0.5, 
            barwidth = unit(10, units = "char"), 
            direction = "horizontal", 
            order = 2
          ),
          size = guide_legend(
            title.position = "top", 
            title.hjust = 0.5, 
            direction = "horizontal", 
            order = 1
          )
        )
    })
    
    caption_plots <- sapply(enrich_sets, FUN = function(.enrich) {
      sprintf(
        fmt = "%s (p-value < %s) and %s (all)", 
        format(length(.enrich@gene), digits = 0, big.mark = ","), 
        pvalue_gene_ora,
        format(length(.enrich@universe), digits = 0, big.mark = ",")
      )
    })
    caption_plots <- sprintf(fmt = "<b>%s</b>) %s",
      LETTERS[seq_along(caption_plots)],
      caption_plots
    )
    
    subtitle_plots <- sprintf(fmt = "<b>%s</b>) %s",
      LETTERS[seq_along(enrich_plots)], 
      names(enrich_plots)
    )

    agg_png(
      filename = file.path(output_directory, "gene_set_enrichment.png"),
      width = 16, height = 16, units = "cm", res = 300, scaling = 0.60
    )
      print(
        wrap_plots(c(enrich_plots, list(guide_area())), ncol = 2, guides = "collect") +
          plot_annotation(
            title = "Over-representation Test",
            subtitle = sprintf("With %s and %s.",
              paste(subtitle_plots[-length(subtitle_plots)], collapse = ", "),
              subtitle_plots[length(subtitle_plots)]
            ),
            caption = sprintf("Mapped genes for %s and %s.",
              paste(caption_plots[-length(caption_plots)], collapse = ", "),
              caption_plots[length(caption_plots)]
            ),
            tag_levels = "A",
            theme = theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
              theme(
                plot.title.position = "plot",
                plot.caption.position = "plot",
                plot.title = element_markdown(),
                plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
                plot.caption = element_markdown(face = "italic", size = rel(0.75)),
                axis.title.x = element_markdown(),
                axis.text.x = element_markdown(),
                axis.title.y = element_markdown(),
                axis.text.y = element_markdown(),
                panel.grid.minor = element_blank()
              )
          )
      )
    invisible(dev.off())
  }
))

## Gene-set Enrichment -----------------------------------------------------------------------------
invisible(lapply(
  X = results_file,
  FUN = function(.file) {
    results <- fread(.file)[order(pvalue)]
    enrich_sets <- list(
      "Reactome" = {
        genes_list <- results[
          pvalue < pvalue_gene & !is.na(entrezgene_id) & entrezgene_id != ""
        ][
          !duplicated(entrezgene_id)
        ][
          i = order(-log2FoldChange),
          j = setNames(log2FoldChange, entrezgene_id)
        ]
        gsePathway(
          geneList = genes_list, 
          organism = organism[["reactome"]],
          pvalueCutoff = fdr_term, 
          pAdjustMethod = "BH"
        )
      },
      "Gene Ontology Biological Process" = {
        genes_list <- results[
          pvalue < pvalue_gene & !is.na(ensembl_gene_id) & ensembl_gene_id != ""
        ][
          !duplicated(ensembl_gene_id)
        ][
          i = order(-log2FoldChange),
          j = setNames(log2FoldChange, ensembl_gene_id)
        ]
        gseGO(
          geneList = genes_list,
          OrgDb = get(organism[["go"]]),
          keyType = "ENSEMBL",
          ont = "BP",
          pvalueCutoff = fdr_term,
          pAdjustMethod = "BH"
        )
      },
      "Gene Ontology Cellular Component" = {
        genes_list <- results[
          pvalue < pvalue_gene & !is.na(ensembl_gene_id) & ensembl_gene_id != ""
        ][
          !duplicated(ensembl_gene_id)
        ][
          i = order(-log2FoldChange),
          j = setNames(log2FoldChange, ensembl_gene_id)
        ]
        gseGO(
          geneList = genes_list,
          OrgDb = get(organism[["go"]]),
          keyType = "ENSEMBL",
          ont = "CC",
          pvalueCutoff = fdr_term,
          pAdjustMethod = "BH"
        )
      },
      "Gene Ontology Molecular Function" = {
        genes_list <- results[
          pvalue < pvalue_gene & !is.na(ensembl_gene_id) & ensembl_gene_id != ""
        ][
          !duplicated(ensembl_gene_id)
        ][
          i = order(-log2FoldChange),
          j = setNames(log2FoldChange, ensembl_gene_id)
        ]
        gseGO(
          geneList = genes_list,
          OrgDb = get(organism[["go"]]),
          keyType = "ENSEMBL",
          ont = "MF",
          pvalueCutoff = fdr_term,
          pAdjustMethod = "BH"
        )
      },
      "KEGG" = {
        genes_list <- results[
          pvalue < pvalue_gene & !is.na(uniprotswissprot) & uniprotswissprot != ""
        ][
          !duplicated(uniprotswissprot)
        ][
          i = order(-log2FoldChange),
          j = setNames(log2FoldChange, uniprotswissprot)
        ]
        gseKEGG(
          geneList = genes_list,
          organism = organism[["kegg"]], 
          keyType = "uniprot", 
          pvalueCutoff = fdr_term, 
          pAdjustMethod = "BH"
        )
      }
    )

    write_xlsx(
      x = setNames(lapply(enrich_sets, FUN = function(.enrich) {
        if (is.null(.enrich) || nrow(.enrich@result) == 0) return(data.frame())
        merge(
          x = setDT(.enrich@result), 
          y = setnames(as.data.table(
            x = sapply(.enrich@geneSets, function(.l) {
              paste(sort(intersect(.l, names(.enrich@geneList))), collapse = "/")
            }), 
            keep.rownames = TRUE
          ), c("ID", "genes_set")),
          by = "ID"
        )[
          j = peripheral_enrichment := {
            peripheral_set <- setdiff(
              unlist(tstrsplit(genes_set, "/"), recursive = TRUE),
              unlist(tstrsplit(core_enrichment, "/"), recursive = TRUE)
            )
            if (length(peripheral_set) == 0) {
              NA_character_
            } else  {
              paste(peripheral_set, collapse = "/")
            }
          },
          by = "ID"
        ]
      }), gsub("Gene Ontology", "GO", names(enrich_sets))), 
      path = file.path(output_directory, "gene_set_enrichment.xlsx")
    )
    
    write_xlsx(
      x = setNames(lapply(enrich_sets, FUN = function(.enrich) {
        if (!is.null(.enrich)) return(data.frame())
        merge(
          x = results,
          y = rbindlist(
            mapply(
              FUN = data.table, 
              entrezgene_id = .enrich@geneSets, 
              gsid = names(.enrich@geneSets), 
              SIMPLIFY = FALSE
            )
          ),
          by = "entrezgene_id"
        )
      }), gsub("Gene Ontology", "GO", names(enrich_sets))), 
      path = file.path(output_directory, "gene_set_list.xlsx")
    )
    
    enrich_sets <- enrich_sets[!sapply(enrich_sets, is.null)]
    
    if (
      all(
        sapply(
          X = enrich_sets, 
          FUN = function(.enrich) nrow(setDT(.enrich@result)[p.adjust < fdr_pathway_gsea])
        ) == 0
      )
    ) {
      return(NULL)
    }
    
    plot_scales <- rbindlist(lapply(enrich_sets, FUN = function(.enrich) {
      setDT(.enrich@result)[
        p.adjust < fdr_pathway_gsea
      ][
        j = list(
          p.adjust = max(round(-log10(c(fdr_pathway_gsea, p.adjust)), digits = 0)),
          setSize = max(c(setSize, 0))
        )
      ]
    }))[j = lapply(.SD, max)] + 1
    
    enrich_plots <- lapply(enrich_sets, FUN = function(.enrich) {
      .dt <- setDT(.enrich@result)[
        p.adjust < fdr_pathway_gsea
      ][
        order(enrichmentScore)
      ][
        j = Description := factor(Description, levels = unique(Description))
      ]
      
      if (nrow(.dt) == 0) {
        return(
          ggplot() +
            annotate(geom = "text", x = 1, y = 1, label = paste0("No enriched terms found\n(FDR < ", fdr_pathway_gsea, ").")) +
            theme(
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              panel.grid = element_blank()
            )
        )
      }
      
      is_too_big <- nrow(.dt) > 25
      
      if (is_too_big) .dt <- first(.dt[order(pvalue)], 25)
      
      ggplot(data = .dt) +
        aes(x = enrichmentScore, y = Description, colour = -log10(p.adjust), size = setSize) +
        geom_point() +
        scale_y_discrete(
          labels = function(.x) {
            .x <- ifelse(
              nchar(.x) > 50,
              gsub("(.{0,50}).*", "\\1 [...]", .x),
              .x
            )
            gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", .x, perl = TRUE)
          }
        ) +
        scale_colour_continuous(limits = c(-log10(fdr_pathway_gsea), plot_scales[["p.adjust"]])) +
        scale_size_continuous(limits = c(1, plot_scales[["setSize"]]), range = c(0.1, 3)) +
        labs(
          x = if (is_too_big) "Enrichment Score (Top 25)" else "Enrichment Score", 
          y = NULL, 
          colour = "-log<sub>10</sub>(FDR)<br>(Benjamini-Hochberg)",
          size = "Number of Genes"
        ) +
        theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
        theme(
          plot.title.position = "plot",
          plot.caption.position = "plot",
          plot.title = element_markdown(),
          plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
          plot.caption = element_markdown(face = "italic", size = rel(0.75)),
          axis.title.x = element_markdown(),
          axis.text.x = element_markdown(),
          axis.title.y = element_markdown(),
          axis.text.y = element_markdown(size = rel(0.75)),
          panel.grid.minor = element_blank(),
          egend.title = element_markdown()
        ) +
        guides(
          colour = guide_colourbar(
            title.position = "top", 
            title.hjust = 0.5, 
            barwidth = unit(10, units = "char"), 
            direction = "horizontal", 
            order = 2
          ),
          size = guide_legend(
            title.position = "top", 
            title.hjust = 0.5, 
            direction = "horizontal", 
            order = 1
          )
        )
    })
    
    caption_plots <- sapply(enrich_sets, FUN = function(.enrich) {
      sprintf(
        fmt = "%s (p-value < %s)", 
        format(length(.enrich@geneList), digits = 0, big.mark = ","), 
        pvalue_gene_gsea
      )
    })
    caption_plots <- sprintf(fmt = "<b>%s</b>) %s",
      LETTERS[seq_along(caption_plots)],
      caption_plots
    )
    
    subtitle_plots <- sprintf(fmt = "<b>%s</b>) %s",
      LETTERS[seq_along(enrich_plots)], 
      names(enrich_plots)
    )
    
    agg_png(
      filename = file.path(output_directory, "gene_set_enrichment.png"),
      width = 16, height = 16, units = "cm", res = 300, scaling = 0.60
    )
      print(
        wrap_plots(c(enrich_plots, list(guide_area())), ncol = 2, guides = "collect") +
          plot_annotation(
            title = "Gene Set Enrichment Analysis",
            subtitle = sprintf("With %s and %s.",
              paste(subtitle_plots[-length(subtitle_plots)], collapse = ", "),
              subtitle_plots[length(subtitle_plots)]
            ),
            caption = sprintf("Mapped genes for %s and %s.",
              paste(caption_plots[-length(caption_plots)], collapse = ", "),
              caption_plots[length(caption_plots)]
            ),
            tag_levels = "A",
            theme = theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
              theme(
                plot.title.position = "plot",
                plot.caption.position = "plot",
                plot.title = element_markdown(),
                plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
                plot.caption = element_markdown(face = "italic", size = rel(0.75)),
                axis.title.x = element_markdown(),
                axis.text.x = element_markdown(),
                axis.title.y = element_markdown(),
                axis.text.y = element_markdown(),
                panel.grid.minor = element_blank()
              )
          )
      )
    invisible(dev.off())
  }
))


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
message(timestamp(quiet = TRUE))
