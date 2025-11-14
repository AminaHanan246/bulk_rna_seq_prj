
# Install packaes if missing; else loads
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# PCA plots
plot_PCA <- function (data,filepath) {
  
  install_if_missing("ggplot2")
  install_if_missing("ggrepel")
  install_if_missing("DESeq2")
  
  png(filename = filepath, 
      width = 2000, height = 2000, res = 300)  # adjust width/height as needed
  
  pcaData <- plotPCA(data,  intgroup = c("condition"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
         title = "PCA Plot colored by condition") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
  print(pca_plot)
  dev.off()
  return(pca_plot)
}


#Disctance plots
plotDists <- function (vsd.obj) {
  
  install_if_missing("pheatmap")
  install_if_missing("RColorBrewer")
  
  png(filename = "results/sampleheatmap.png", width = 1000, height = 900, res = 300)
  
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( vsd.obj$condition )
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(55)
  plot<- pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors,
                     fontsize_col = 4,
                     fontsize_row = 4,
                     fontsize_legend = 4,
                     fontsize = 4)
  print(plot)
  dev.off()
  return(plot)
}

# Variable gene Heatmap
variable_gene_heatmap <- function (vsd.obj, num_genes = 500, annotation, title = "") {
  
  install_if_missing("pheatmap")
  install_if_missing("RColorBrewer")
  install_if_missing("matrixStats")
  
  png(filename = "results/variable_gene_heatmap.png", 
      width = 2000, height = 1500, res = 300)  # adjust width/height as needed
  
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  stabilized_counts <- assay(vsd.obj)
  row_variances <- rowVars(stabilized_counts)
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing=T)[1:num_genes],]
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm=T)
  gene_names <- annotation$Gene.name[match(rownames(top_variable_genes), annotation$Gene.stable.ID)]
  rownames(top_variable_genes) <- gene_names
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL
  heatmap<- pheatmap::pheatmap(top_variable_genes, color = mr, annotation_col = coldata, fontsize_col = 4, fontsize_row = 250/num_genes, border_color = NA, main = title)
  print(heatmap)
  dev.off()
  return(heatmap)
  }

# Gene expression across condition
gene_plot_counts <- function (dds, gene, normalization = "DESeq2"){
  
  install_if_missing("tibble")
  install_if_missing("ggplot2")
  install_if_missing("edgeR")
  
  if (normalization == "cpm") {
    normalized_data <- cpm(counts(dds, normalized = F)) 
  } else if (normalization == "DESeq2")
    normalized_data <- counts(dds, normalized = T) 
  condition <- dds@colData$condition
  if (is.numeric(gene)) { 
    if (gene%%1==0 )
      ensembl_id <- rownames(normalized_data)[gene]
    else
      stop("Invalid index supplied.")
  } else if (gene %in% annotation$Genesymbol){ 
    ensembl_id <- annotation$Geneid[which(annotation$Genesymbol == gene)]
  } else if (gene %in% annotation$Geneid){
    ensembl_id <- gene
  } else {
    stop("Gene not found. Check spelling.")
  }
  
  expression <- normalized_data[ensembl_id,]
  gene_name <- annotation$Genesymbol[which(annotation$Geneid == ensembl_id)]
  gene_tib <- tibble(condition = condition, expression = expression)
  gene_plot <- ggplot(gene_tib, aes(x = condition, y = expression))+
    geom_boxplot(outlier.size = NULL)+
    geom_point()+
    labs (title = paste0("Expression of ", gene_name, " - ", ensembl_id), x = "group", y = paste0("Normalized expression (", normalization , ")"))+
    theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11))
  ggsave(filename =paste0("results/",gene,"_cond.png") , plot = gene_plot,bg = "white", width = 8, height = 6, dpi = 300)
  return(gene_plot)
}
