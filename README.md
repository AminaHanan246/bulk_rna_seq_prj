# Bulk RNA-Seq Analysis Project with GSEA

This repository contains a complete workflow for bulk RNA-seq analysis, from raw sequence data to differential gene expression and pathway enrichment analysis. The project demonstrates RNA-seq preprocessing, quality control, alignment, quantification, differential expression, and Gene Set Enrichment Analysis (GSEA) using R, Python, and standard bioinformatics tools.

---

Project Overview

The goal of this project is to analyze bulk RNA-seq data to identify differentially expressed genes (DEGs) and enriched biological pathways. The workflow includes:

1. Data Acquisition: Download raw RNA-seq data (SRA format) and convert to FASTQ.
2. Quality Control: Assess read quality using FastQC and MultiQC.
3. Preprocessing: (Optional) trimming of adapters and low-quality reads using Trimmomatic.
4. Alignment: Map reads to reference genome using HISAT2 to human reference genome GRCh38_p14 .
5. Quantification: Count reads per gene using featureCounts.
6. Differential Expression Analysis: Perform DE analysis with DESeq2.
7. Gene Set Enrichment Analysis (GSEA): Identify enriched pathways from DE genes.
8. Visualization: Generate PCA plots, heatmaps, volcano plots, and GSEA enrichment plots.

---

Technologies & Tools

* Languages: R, Python, Bash
* Bioinformatics Tools: FastQC, MultiQC, Trimmomatic, HISAT2, featureCounts, DESeq2, GSEA
* Visualization: ggplot2, pheatmap, waterfallplots
* Data Analysis: Differential expression, PCA, clustering, pathway enrichment
* Version Control: Git & GitHub

---

Folder Structure

bulk\_rna\_seq\_prj/
|
+-- data/               # Processed RNA-seq data and metadata
+-- scripts/            # Analysis scripts (R, Python, Bash)
+-- results/            # DE results, GSEA results, plots
+-- README.txt          # Project documentation
+-- environment.yml     # dependencies for reproducibility

---

Setup & Installation

1. Clone the repository:

```
git clone https://github.com/AminaHanan246/bulk_rna_seq_prj.git
cd bulk_rna_seq_prj
```

2. Install dependencies:

* R packages: DESeq2, ggplot2, pheatmap, clusterProfiler (for GSEA)
* Python packages: pandas
* Bioinformatics tools: FastQC, Trimmomatic, HISAT2, featureCounts

---

Usage

1. Downloading SRA files:

```
python3 scripts/fastq_download.py
```

2. Trimming, file handling and HISAT2 alignment:

```
python3 scripts/hisat_align.py
```

3. Quantification of reads:

```
python3 scripts/feature_counts.py
```

3. Differential Expression Analysis & GSEA:

```
Rscript scripts/bulk_rna_seq_GSE106305.Rmd
```

5. Visualization:

* PCA plots, heatmaps, volcano plots, and GSEA enrichment plots are saved in results/.

---

Key Learnings

* Familiarity with NGS pipelines: QC, trimming, alignment, quantification, differential expression
* Experience with GSEA for pathway and functional enrichment analysis
* Hands-on experience with R for statistical analysis and visualization
* Ability to handle real RNA-seq datasets and interpret biological meaning

---
