# Bulk RNA-Seq Analysis Project with GSEA

This repository contains a complete workflow for bulk RNA-seq analysis, from raw sequence data to differential gene expression and pathway enrichment analysis. The project demonstrates RNA-seq preprocessing, quality control, alignment, quantification, differential expression, and Gene Set Enrichment Analysis (GSEA) using R, Python, and standard bioinformatics tools.

---

## Project Overview

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

## Folder Structure
=======

bulk\_rna\_seq\_prj/
|
+-- data/               # Processed RNA-seq data and metadata
+-- scripts/            # Analysis scripts (R, Python, Bash)
+-- results/            # DE results, GSEA results, plots
+-- README.txt          # Project documentation
+-- environment.yml     # dependencies for reproducibility

---
## Setup & Installation
=======
**Programs required:** it is recommended that the user has anaconda installed, through which all required programs can be installed. Assuming that anaconda is available, all the required programs can be installed using the following:
```bash
#Install the required programs using anaconda
conda create -n preprocess python=3.7

conda install -n preprocess -c bioconda fastqc
conda install -n preprocess -c trimmomatic
conda install -n RNA-seq -c bioconda multiqc
conda install -n preprocess -c hisat2
conda install -n preprocess -c samtools
```
The NCBI’s SRA toolkit is installed, and the path is added to /.bashrc configuration file:
```bash
#Install SRA toolkit
sudo apt install sra-toolkit

#Adding PATH to configuration file
nano ~/.bashrc
export PATH="/home/amina/sratoolkit/bin:$PATH"
source ~/.bashrc
```
The conda environment needs to be activated before running the pipeline using command:
```bash
conda activate preprocess
```
---
## Pre-processing Pipeline
=======

The dataset produced as part of the research [Guo et al. Nature Communications 2019](https://www.ncbi.nlm.nih.gov/pubmed/30655535). The raw sequencing data are obtained via Gene Expression Omnibus (GEO) using the accession number provided in the publication : GSE106305. The corresponding webpage contains information regarding the sequencing data for each sample in the study. 
The objective was to find the differentially expressed under conditions: normoxia, normal oxygen condition of typically 21% O2, and hypoxia, low oxygen condition of typically 1-5% O2, for cell lines LNCAP and PC3. For this the control samples (Empty\_Vector for LNCaP and siCtrl for PC3) where selected. The samples to be downloaded and associated SRA accession number are in the table below:
| Sample Name                                   | GSM Identifier | SRA Identifier (SRX) | SRA Runs (SRR, download these)                     |
|-----------------------------------------------|----------------|----------------------|----------------------------------------------------|
| LNCaP\_RNA-Seq\_Empty\_Vector\_Normoxia\_rep1 | GSM3145509     | SRX4096735           | SRR7179504, SRR7179505, SRR7179506, and SRR7179507 |
| LNCaP\_RNA-Seq\_Empty\_Vector\_Normoxia\_rep2 | GSM3145510     | SRX4096736           | SRR7179508, SRR7179509, SRR7179510, and SRR7179511 |
| LNCaP\_RNA-Seq\_Empty\_Vector\_Hypoxia\_rep1  | GSM3145513     | SRX4096739           | SRR7179520, SRR7179521, SRR7179522, and SRR7179523 |
| LNCaP\_RNA-Seq\_Empty\_Vector\_Hypoxia\_rep2  | GSM3145514     | SRX4096740           | SRR7179524, SRR7179525, SRR7179526, and SRR7179527 |
| PC3\_RNA-Seq\_siCtrl\_Normoxia\_rep1          | GSM3145517     | SRX4096743           | SRR7179536                                         |
| PC3\_RNA-Seq\_siCtrl\_Normoxia\_rep2          | GSM3145518     | SRX4096744           | SRR7179537                                         |
| PC3\_RNA-Seq\_siCtrl\_Hypoxia\_rep1           | GSM3145521     | SRX4096747           | SRR7179540                                         |
| PC3\_RNA-Seq\_siCtrl\_Hypoxia\_rep2           | GSM3145522     | SRX4096748           | SRR7179541                                         |
### Download FASTQ files using SRA tools
The sequencing data stored in NCBI database in the form SRA files, which can be downloaded using NCBI’s SRA toolkit. The `prefetch` command downloads the SRA files based on the specified SRA accession number. 
```bash
prefetch SRR7179504

2025-08-27T05:28:06 prefetch.3.2.1: 1) Resolving 'SRR7179504'...
2025-08-27T05:28:10 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-08-27T05:28:12 prefetch.3.2.1: 1) Downloading 'SRR7179504'...
2025-08-27T05:28:12 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-08-27T05:28:12 prefetch.3.2.1:  Downloading via HTTPS...
2025-08-27T05:29:00 prefetch.3.2.1:  HTTPS download succeed
2025-08-27T05:29:04 prefetch.3.2.1:  'SRR7179504' is valid: 439677804 bytes were streamed from 439667257
2025-08-27T05:29:04 prefetch.3.2.1: 1) 'SRR7179504' was downloaded successfully
2025-08-27T05:29:04 prefetch.3.2.1: 1) Resolving 'SRR7179504's dependencies...
2025-08-27T05:29:04 prefetch.3.2.1: 'SRR7179504' has 0 unresolved dependencies
```
The downloaded SRA file “SRR7019504” is then read into FASTQ file. For this the command `fastq-dump`:
```bash
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip /mnt/d/BI_prj/bulkrnaseq_proj/normoxia_vs_hypoxia/SRR7179504/SRR7179504.sra
Rejected 13548432 READS because of filtering out non-biological READS
Read 13548432 spots for /mnt/d/BI_prj/bulkrnaseq_proj/normoxia_vs_hypoxia/SRR7179504/SRR7179504.sra
Written 13548432 spots for /mnt/d/BI_prj/bulkrnaseq_proj/normoxia_vs_hypoxia/SRR7179504/SRR7179504.sra
```
A compressed file `SRR7179504_pass.fastq.gz` is created in the subdirectory called `fastq`. 
Since multiple SRA files are to be downloaded, a python script is written to automate the process. The code is provided in [`scripts/fastq_download.py`](scripts/fastq_download.py) and is as follows:
```python
import subprocess
import time

sra_numbers = [
    "SRR7179504", "SRR7179505", "SRR7179506", "SRR7179507",
    "SRR7179508", "SRR7179509", "SRR7179510", "SRR7179511",
    "SRR7179520", "SRR7179521", "SRR7179522", "SRR7179523",
    "SRR7179524", "SRR7179525", "SRR7179526", "SRR7179527",
    "SRR7179536", "SRR7179537", "SRR7179540", "SRR7179541"
]

for sra_id in sra_numbers:
    print("\n=== Downloading:", sra_id, "===")
    prefetch_cmd = f"prefetch {sra_id}"
    print("Command:", prefetch_cmd)

    start_time = time.time()
    subprocess.call(prefetch_cmd, shell=True)
    end_time = time.time()

    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ Download time for {sra_id}: {elapsed_min:.2f} minutes")

for sra_id in sra_numbers:
    sra_path = f'"/mnt/d/BI_prj/bulkrnaseq_proj/normoxia_vs_hypoxia/data_analysis/{sra_id}/{sra_id}.sra"'
    print("\n=== Generating FASTQ for:", sra_id, "===")
    fastq_dump_cmd = (
        f"fastq-dump --outdir fastq --gzip --skip-technical "
        f"--readids --read-filter pass --dumpbase --split-3 --clip {sra_path}"
    )
    print("Command:", fastq_dump_cmd)

    start_time = time.time()
    subprocess.call(fastq_dump_cmd, shell=True)
    end_time = time.time()

    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ FASTQ generation time for {sra_id}: {elapsed_min:.2f} minutes")
```
All the SRA files are now downloaded in subdirectory `fastq` wth path to sra file being `/SRR*/SRR*.sra`
### Pre-alignment QC
The raw sequence data is assessed for quality. FastQC reports are generated for samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination with the results output to the subdirectory`/fastq_results` using the following command:
```bash
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
```
The fastqc reports can be combined into one summary report using `Mulitqc` with the following command:
```bash
multiqc fastqc_results/ -o multiqc_report/
```
### Trimming(optional)
Trimming is pre-alignment step to remove adapter sequences and low-quality bases. The step based on the fastqc report generated earlier which show the quality scores. The following command is used to trim SRR7079504:
```bash
trimmomatic SE -threads 4 SRR7179504_pass.fastq.gz SRR7179504_trimmed.fastq.gz TRAILING:10 -phred33
```
The quality of the read is checked again after the trimming.
Since the reads had quality scores above 10 and adapter sequences removed during FASTQ file conversion, this step is optional (lenient filtering is done since the reads are used for quantification).

>>>>>>> 35ed87f (Updating README and scripts)

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

## Usage

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

## Key Learnings

* Familiarity with NGS pipelines: QC, trimming, alignment, quantification, differential expression
* Experience with GSEA for pathway and functional enrichment analysis
* Hands-on experience with R for statistical analysis and visualization
* Ability to handle real RNA-seq datasets and interpret biological meaning

---
