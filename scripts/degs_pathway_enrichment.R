

# Install packaes if missing; else loads
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Filtering data based on cell type
filter_celltype <- function (dds,condition,ref = "Normoxia",against = "Hypoxia"){
  # Condition labels
  ref_label <- paste0(toupper(cell_type),"_",ref)
  contrast_label <- paste0(toupper(cell_type),"_",against)
  
  # filter dds based on cell type
  dds_cell_filtered <- dds[, grepl(paste0(toupper(cell_type)), colnames(dds))]
  dds_cell_filtered
  dds_cell_filtered$condition <- droplevels(dds_cell_filtered$condition)#removing unrelated levels
  
  #setting up reference
  dds_cell_filtered$condition <- relevel(dds_cell_filtered$condition, ref = ref_label)
  dds_cell_filtered <- DESeq(dds_cell_filtered)
  
  # Extract differential expression results for contrast: Hypoxia vs Normoxia
  res_celltype <- results(dds_cell_filtered, contrast = c("condition", contrast_label, ref_label))
  res_celltype
  summary(res_celltype)
  
  res_celltypeOrdered <- res_celltype[order(res_celltype$padj), ] #order with padj values
  sum(res_celltypeOrdered$padj < 0.05, na.rm = TRUE)
  head(res_celltypeOrdered)
  summary(res_celltypeOrdered)
  write.csv(as.data.frame(res_celltypeOrdered), file = paste0("data/DEGs_",cell_type,".csv"))
  
  return(list(
    dds = dds_cell_filtered,
    results = res_celltype,
    ordered_results = res_celltypeOrdered
  )) 
}

# Volcano plot - gene regulation
gene_regulation_plot <- function(ordered_results){
  install_if_missing("ggplot2")
  
  res_df <- as.data.frame(ordered_results) 
  res_df <- na.omit(res_df)
  res_df$gene <- rownames(res_df) #gene names added as column
  
  #Gene regulation categorised
  res_df$regulation <- "Not Significant" #defaulted 
  res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
  res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"
  
  #Volcano plot
  gene_plot <- ggplot(res_df, aes(x = log2FoldChange, 
                                  y = -log10(padj), 
                                  color = regulation)) + #categorise based on regulation
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "#FEA405", 
                                  "Downregulated" = "purple", 
                                  "Not Significant" = "gray")) +
    #Threshold line at padj = 0.05
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    #placement of text
    annotate("text", x = min(res_df$log2FoldChange), y = -log2(0.05) + 0.5,
             label = "padj = 0.05", hjust = 0, size = 3) +
    theme_minimal() +
    #labels
    labs(title = paste0("Volcano Plot -",toupper(cell_type)), 
         x = "Log2 Fold Change", 
         y = "-Log10 Adjusted P-Value") +
    theme(plot.title = element_text(hjust = 0.5))
  v_plot <- paste0("results/vp_",cell_type,".png")
  ggsave(v_plot, plot = gene_plot,bg = "white", width = 8, height = 6, dpi = 300)
  
  return(gene_plot)
}

# Gene Set Enrichment Analysis (GSEA)
convert_to_entrez <- function(results){
  install_if_missing("clusterProfiler")
  install_if_missing("org.Hs.eg.db")
  install_if_missing("dplyr")
  install_if_missing("stats")
  
  # Convert DESeqResults to data.frame
  results <- as.data.frame(results)
  results$ENSEMBL <- rownames(results)
  
  #convert ENSEMBL ids to entrez ids for Reactome
  ncbi_list <- clusterProfiler::bitr(
    geneID = rownames(results),        # use Ensembl IDs from row names
    fromType = "ENSEMBL",          
    toType = "ENTREZID", 
    OrgDb = org.Hs.eg.db
  )
  #Column with Esemble ids
  results$ENSEMBL <- rownames(results)
  head(results)
  
  #Merge ncbi_list with valid entrez id and remove duplicates 
  res_mapped <- results %>%
    left_join(ncbi_list, by = "ENSEMBL") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  #Genes ranked based on log2FoldChange
  ngenes <- res_mapped$log2FoldChange
  names(ngenes) <- res_mapped$ENTREZID
  ngenes <- sort(ngenes, decreasing = TRUE)
  
  return(list(
    res_mapped = res_mapped,
    ngenes = ngenes
  ))
}

#GSEA
run_gsea <- function(ngenes) {  # Fixed: created proper function
  install_if_missing("ReactomePA")
  gsea <- gsePathway(
    ngenes,
    organism = "human",
    verbose = FALSE
  )
  return(gsea)
}


top_20_pathways <- function(gsea){
  install_if_missing("dplyr")
  install_if_missing("forcats")
  
  pathways <- gsea@result 
  # Sort pathways by adjusted p-value (FDR) to prioritize statistically significant ones
  pathways <- pathways[order(pathways$p.adjust), ]
  
  #Reorder by absolute NES to highlight strongest up/down regulated pathways
  top_pathways <- pathways[order(abs(pathways$NES), decreasing = TRUE), ]  # Sort by NES
  
  top20 <- top_pathways[1:20, ] %>%
    mutate(Description = fct_reorder(Description, NES)) # Reorder factor for y-axis based on NES
  
  write.csv(top20, paste0("data/",cell_type,"_top20_pathways.csv"), row.names = FALSE)
  return(top20)
}

# Pathway enrichment based on overall expression profile

pathway_enrich_overall <- function(top20){ 
  install_if_missing("ggplot2")
  
  enrich_overall <- ggplot(top20, aes(x = NES,                 
                                      y = Description,
                                      color = p.adjust,
                                      size = setSize)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#0072B2", high = "#D55E00", name = "FDR (p.adjust)") +
    scale_size(range = c(3, 10), name = "Gene Set Size") +
    labs(
      title = paste0("Top 20 Overall Enriched Pathways in ",cell_type),
      subtitle = "Gene Set Enrichment Analysis (GSEA)",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      caption = "Data source: clusterProfiler::gsePathway"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 13),
      legend.position = "right"
    )
  ggsave(paste0("results/enrichment_overall_",cell_type,".png"), 
         plot = enrich_overall,
         bg = "white",
         width = 15, height = 6, dpi = 300)
  
  return(enrich_overall)
}

# Pathway enrichment by significant DEGs
pathway_enrich_sigDEGS <- function(res_mapped, organism = "human"){
  install_if_missing("ReactomePA")
  
  # Filter DEGs: padj < 0.1 and |log2FC| > 0.5, then extract ENTREZ IDs
  sig_genes <- res_mapped %>%
    filter(padj < 0.1, abs(log2FoldChange) > 0.5) %>%
    pull(ENTREZID)
  
  # Run Reactome pathway enrichment for significant genes
  enr <- enrichPathway(gene = sig_genes, organism = organism, pvalueCutoff = 0.1)
  
  reactome_plot <- dotplot(enr, showCategory=20)+
    labs(
      title = paste0("Top 20 DEGs Enriched Pathways in ",toupper(cell_type)))
  ggsave(paste0("results/react_sig_genes_",cell_type,".png"), 
         plot = reactome_plot, 
         bg = "white",
         width = 8, height = 10, dpi = 300)
  
  return(reactome_plot)
}


# GSEA of Hallmark Programs 

# Map Ensembl IDs to gene symbols using org.Hs.eg.db
ranked_gene_list <- function(res_celltype){
  
  res_celltype$ENSEMBL <- rownames(res_celltype)
  
  symbol_map <- bitr(res_celltype$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  
  # Merge and filter
  rank_df <- merge(res_celltype, symbol_map, by.x = "ENSEMBL", by.y = "ENSEMBL")
  rank_df <- rank_df[, c("SYMBOL", "log2FoldChange")]
  colnames(rank_df) <- c("Gene.name", "log2FoldChange")
  
  # Save as .rnk
  head(rank_df)
  write.table(rank_df, file = "data/lncaprank.rnk", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  return(rank_df)
}


# Prepare and validate Ranked Gene list
prepare_ranked_list <- function(ranked_list) {
  #check if nammed numeric vector
  if (is.vector(ranked_list) && !is.list(ranked_list)) {
    return(ranked_list)  # Return as-is if already processed
  }
  #check for requiered columns
  if (!is.data.frame(ranked_list)) {
    stop("Input 'ranked_list' must be a data frame with 'Gene.name' and 'log2FoldChange' columns.")
  }
  #Check duplicate genes and avg. fold change
  if (sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(. ~ Gene.name, data = ranked_list, FUN = mean)
    ranked_list <- ranked_list[order(ranked_list$log2FoldChange, decreasing = TRUE), ]
  }
  
  ranked_list <- na.omit(ranked_list)
  ranked_list <- tibble::deframe(ranked_list[, c("Gene.name", "log2FoldChange")])
  
  return(ranked_list)
}

run_fgsea <- function(hallmark_pathway, ranked_list) {
  install_if_missing("fgsea")
  install_if_missing("stats")
  
  fgsea_results <- fgsea(pathways = hallmark_pathway,
                         stats = ranked_list,
                         minSize = 15,
                         maxSize = 500,
                         nperm = 1000)
  
  fgsea_results_ordered <- fgsea_results[order(-fgsea_results$NES), ]
  
  return(list(
    results = fgsea_results,
    ordered_results = fgsea_results_ordered
  ))
}

# Hallmark – Cell Hypoxia Response
# Waterfall plot for Hallmark pathways
waterfall_plot <- function(fgsea_results, graph_title) {
  install_if_missing("stringr")
  install_if_missing("ggplot2")
  install_if_missing("dplyr")
  
  if (is.null(fgsea_results)) {
    stop("fgsea_results is NULL. Cannot generate plot.")
  }
  
  plot_obj <- fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_", 2)[, 2]) %>%
    ggplot(aes(reorder(short_name, NES), NES)) +
    geom_bar(stat = "identity", aes(fill = padj < 0.05)) +
    coord_flip() +
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7), 
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("results/hallmark_enrich_", cell_type, ".png"),
         plot = plot_obj,
         bg = "white",
         width = 10, height = 10, dpi = 300)
  return(plot_obj)
}
