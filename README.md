# Bulk RNA-Seq Analysis Project with GSEA

This repository contains a complete workflow for bulk RNA-seq analysis, from raw sequence data to differential gene expression and pathway enrichment analysis. The project demonstrates RNA-seq preprocessing, quality control, alignment, quantification, differential expression, and Gene Set Enrichment Analysis (GSEA) using R, Python, and standard bioinformatics tools.

---

## Project Overview

The goal of this project is to analyze bulk RNA-seq data to identify differentially expressed genes (DEGs) and enriched biological pathways. The workflow includes:

1. Data Acquisition: Download raw RNA-seq data (SRA format) and convert to FASTQ.
2. Quality Control: Assess read quality using FastQC and MultiQC.
3. Pre-processing: (Optional) trimming of adapters and low-quality reads using Trimmomatic.
4. Alignment: Map reads to reference genome using HISAT2 to human reference genome GRCh38_p14 .
5. Sorting and indexing using SAMtools 
6. Quantification: Count reads per gene using featureCounts.
7. Differential Expression Analysis: Perform DE analysis with DESeq2.
8. Gene Set Enrichment Analysis (GSEA): Identify enriched pathways from DE genes.
9. Visualization: Generate PCA plots, heatmaps, volcano plots, and GSEA enrichment plots.

---

## Folder Structure
```
bulk_rna_seq_prj/
├── data/               # Processed RNA-seq data and metadata
├── scripts/            # Analysis scripts (R, Python, Bash)
├── results/            # DE results, GSEA results, plots
├── README.txt          # Project documentation
└── environment.yml     # Dependencies for reproducibility
```
---
## Setup & Installation

**Programs required:** It is recommended that the user have Anaconda installed, through which all required programs can be installed. Assuming that Anaconda is available, all the required programs can be installed using the following:
```bash
#Install the required programs using anaconda
conda create -n preprocess python=3.7

conda install -n preprocess -c bioconda fastqc
conda install -n preprocess -c trimmomatic
conda install -n RNA-seq -c bioconda multiqc
conda install -n preprocess -c hisat2
conda install -n preprocess -c samtools
```

The NCBI’s SRA toolkit is installed, and the path is added to the /.bashrc configuration file:
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

The dataset was produced as part of the research [Guo et al., Nature Communications 2019](https://www.ncbi.nlm.nih.gov/pubmed/30655535). The raw sequencing data are obtained via Gene Expression Omnibus (GEO) using the accession number provided in the publication: GSE106305. The corresponding webpage contains information regarding the sequencing data for each sample in the study. 
The objective was to find the differentially expressed under conditions of normoxia, normal oxygen condition of typically 21% O2, and hypoxia, low oxygen condition of typically 1-5% O2, for cell lines LNCAP and PC3. For this, the control samples (Empty\_Vector for LNCaP and siCtrl for PC3) were selected. The samples to be downloaded and associated SRA accession numbers are in the table below:
| Sample Name                                   | GSM Identifier | SRA Identifier (SRX) | SRA Runs (SRR, download these)                     |
|-----------------------------------------------|----------------|----------------------|----------------------------------------------------|
| LNCaP\_RNA-Seq\_Empty\_Vector\_Normoxia\_rep1 | GSM3145509     | SRX4096735           | SRR7179504  SRR7179505  SRR7179506  SRR7179507 |
| LNCaP\_RNA-Seq\_Empty\_Vector\_Normoxia\_rep2 | GSM3145510     | SRX4096736           | SRR7179508  SRR7179509  SRR7179510  SRR7179511 |
| LNCaP\_RNA-Seq\_Empty\_Vector\_Hypoxia\_rep1  | GSM3145513     | SRX4096739           | SRR7179520  SRR7179521  SRR7179522  SRR7179523 |
| LNCaP\_RNA-Seq\_Empty\_Vector\_Hypoxia\_rep2  | GSM3145514     | SRX4096740           | SRR7179524  SRR7179525  SRR7179526  SRR7179527 |
| PC3\_RNA-Seq\_siCtrl\_Normoxia\_rep1          | GSM3145517     | SRX4096743           | SRR7179536                                         |
| PC3\_RNA-Seq\_siCtrl\_Normoxia\_rep2          | GSM3145518     | SRX4096744           | SRR7179537                                         |
| PC3\_RNA-Seq\_siCtrl\_Hypoxia\_rep1           | GSM3145521     | SRX4096747           | SRR7179540                                         |
| PC3\_RNA-Seq\_siCtrl\_Hypoxia\_rep2           | GSM3145522     | SRX4096748           | SRR7179541                                         |

Download FASTQ files using SRA tools
------------------------------------
The sequencing data is stored in the NCBI database in the form of SRA files, which can be downloaded using NCBI’s SRA toolkit. The `prefetch` command downloads the SRA files based on the specified SRA accession number. 
```bash
prefetch SRR7179504

2025-08-27T05:28:06 prefetch.3.2.1: 1) Resolving 'SRR7179504'...
2025-08-27T05:28:10 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-08-16T05:28:12 prefetch.3.2.1: 1) Downloading 'SRR7179504'...
2025-08-16T05:28:12 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-08-16T05:28:12 prefetch.3.2.1:  Downloading via HTTPS...
2025-08-16T05:29:00 prefetch.3.2.1:  HTTPS download succeed
2025-08-16T05:29:04 prefetch.3.2.1:  'SRR7179504' is valid: 439677804 bytes were streamed from 439667257
2025-08-16T05:29:04 prefetch.3.2.1: 1) 'SRR7179504' was downloaded successfully
2025-08-16T05:29:04 prefetch.3.2.1: 1) Resolving 'SRR7179504's dependencies...
2025-08-16T05:29:04 prefetch.3.2.1: 'SRR7179504' has 0 unresolved dependencies
```
The downloaded SRA file “SRR7019504” is then read into FASTQ file. For this the command `fastq-dump`:
```bash
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip /mnt/d/BI_prj/bulkrnaseq_proj/normoxia_vs_hypoxia/SRR7179504/SRR7179504.sra

Rejected 13548432 READS because of filtering out non-biological READS
Read 13548432 spots for /mnt/d/BI_prj/bulkrnaseq_proj/normoxia_vs_hypoxia/SRR7179504/SRR7179504.sra
Written 13548432 spots for /mnt/d/BI_prj/bulkrnaseq_proj/normoxia_vs_hypoxia/SRR7179504/SRR7179504.sra
```
`--outdir fastq`           : Specifies the output directory for the FASTQ files  
`--gzip`                   : Compresses the output FASTQ files using gzip  
`--skip-technical`         : Skips technical reads (e.g., control reads or adapters)  
`--readids`                : Includes read identifiers in the FASTQ header  
`--read-filter pass`       : Filters out low-quality reads; keeps only those marked "pass"  
`--dumpbase`               : Outputs base calls (A, T, G, C, N) instead of color space  
`--split-3`                : Splits paired-end reads into separate files (`_1.fastq.gz`, `_2.fastq.gz`)  
`--clip`                   : Removes adapter sequences from reads  
`~/sra/...`                : Path to the input SRA file  

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
    sra_path = f'"../{sra_id}/{sra_id}.sra"'
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

Pre-alignment QC
----------------
The raw sequence data is assessed for quality. The following command is used to generate FastQC reports for samples in order to evaluate sequence quality, GC content, duplication rates, length distribution, K-mer content, and adapter contamination. The results are then recorded in the subdirectory "/fastq_results."
```bash
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
```

The fastqc reports can be combined into one summary report using `Mulitqc` with the following command:
```bash
multiqc fastqc_results/ -o multiqc_report/
```

Trimming(optional)
------------------
Trimming is pre-alignment step to remove adapter sequences and low-quality bases. The step us done based on the FASTqc report generated earlier which show the quality scores. The following command is used to trim SRR7079504:
```bash
trimmomatic SE -threads 4 SRR7179504_pass.fastq.gz SRR7179504_trimmed.fastq.gz TRAILING:10 -phred33
```
After the trimming process, the quality of the read is once again assessed.

Since the reads have short length and adapter sequences removed during FASTQ file conversion, this step is skipped as it may introduce biasness and poor alignment (lenient filtering is done since the reads are used for quantification).

### Concatenating FASTQ files to sample files
The LNCAP samples are associated with four SRA files each and is concatenated into single FASTQ file using the command:
```bash
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz  > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz  > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz  > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz  > LNCAP_Hypoxia_S2.fastq.gz
```

The PC3 sample are associated with only one SRA file and is therefore renamed to the sample names using command:
```bash
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
```

Mapping reads using HISAT2
---------------------------
The reads from FASTQ file are aligned to a reference genome, where the reads are matched based on sequence similarity in the reference genome. This tells us which part of the gene was transcribed for the mRNA, and the number of times a read is mapped to a specific gene indicates whether the gene expression was high or low.

### Concept behind mapping
To perform bulk-RNA sequence analysis,  library is created. The mRNA transcripts from cells are reverse transcribed into cDNA, fragmented and are attached with adapter sequences in both ends and this created the library. The sequencer reads the fragments and stores the sequence in FASTQ files. The FASTQ file consists of 4-line chunks, as shown
```bash
zcat LNCAP_Hypoxia_S1.fastq.gz | head -4

@SRR7179520.1.1 1 length=76    
GTGAANATAGGCCTTAGAGCACTTGANGTGNTAGNGCANGTNGNNCCGGAACGNNNNNNNNAGGTNGNNNGNGTTG
+SRR7179520.1.1 1 length=76
AAAAA#EEEEEEEEEEEEEEEEEEEE#EEE#EEE#EEE#EE#E##EEEEEEEE########EEEE#E###E#EAEA
```
 1. Sequence identifier(starts with @)
 2. Read sequence
 3. Separator line(starts with +)
 4. Quality score of each base (based on ASCII)

### Mapping reads
The pre-built genome index, required for mapping, is downloaded using the `wget` command and is extracted using `tar -xvzf` command:
```bash
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

tar -xvzf grch38_genome.tar.gz
```
The reads from FASTQ files are aligned to genome index using HISAT2, which does splice-aware alignment for the LNCAP_Normoxia_S1:
```bash
hisat2 -p 8 -x GRCh38_index -U LNCAP_Normoxia_S1_R1_001.fastq.gz -S LNCAP_Normoxia_S1.sam
```
> [!NOTE] 
> STAR aligner provides more accurate and sensitive mapping; however, HISAT2 is used here because it uses less memory and is significantly faster

The output BAM file consists of 11 fields for alignment information: 
```bash
samtools view LNCAP_Normoxia_S2.bam | head -n 1

SRR7179504.1361607.1    272     1       14277   0       76M     *0       0       GAAACAGGGCCGCGGGGAGCGGCTGCCCCCACTGCCTAGGGACCAACAGGGGCAGGAGGCAGTCACTGACCCCGAG     //</EEEE/E//</6//EEE//E/AEEEE<//AEEEAAEEEEEEE<AEEEEEAEEEEEEEEE6EEEEEEEEAAAAA     AS:i:-15        ZS:i:-15 XN:i:0  XM:i:5  XO:i:0  XG:i:0  NM:i:5  MD:Z:8T1A2T6T0T54YT:Z:UU NH:i:3
```
 1. QNAME -Read name
 2. FLAG - bitwise flag
 3. RNAME - Chromosome (if read is not aligned '*')
 4. POS - 1-based left-most mapping position
 5. MAPQ - Mapping quality
 6. CIGAR - describes the position of instertion(I),deletion(D) and matches(M) in alignment
 7. RNEXT - name of the pair sequence (in pair-ended sequences) - * unpaired
 8. PNEXT - position of pair 
 9. TLEN - Total span (size of pair reads and distance between them)
 10. SEQ - Read Sequence
 11. QUAL - Phred quality
 and TAG - TAG information (AS - alignment score, NH - number of reported alignments 

Sorting and Indexing BAM Files using SAMtool
--------------------------------------------
After alignment, the BAM files created are unsorted. The samtools sort and index the BAM file based on genomic coordinates, which is essential for downstream analysis. The aim is to achieve fast retrieval of alignments mapped to regions without having to process the entire set of alignments. The file is sorted based on:
1. Chromosome name
2. Start position 

Since multiple files need to be aligned to the genome index, sorted, and indexed using samtool, the process is automated using a Python script. The code is provided in [`scripts/hisat_align.py`](scripts/hisat_align.py):
```python
files = [
    "LNCAP_Normoxia_S1.fastq.gz",
    "LNCAP_Normoxia_S2.fastq.gz",
    "LNCAP_Hypoxia_S1.fastq.gz",
    "LNCAP_Hypoxia_S2.fastq.gz",
    "PC3_Normoxia_S1.fastq.gz",
    "PC3_Normoxia_S2.fastq.gz"
]

logfile = "alignment_log.txt"

for file in files:
    sample_name = os.path.basename(file).replace(".fastq.gz", "")
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"{timestamp} - Processing {sample_name}\n"
    with open(logfile, "a") as log:
        log.write(log_entry)
    
    start_time = time.time()

    hisat2_exec = (
        f"hisat2 -x grch38_genome -U fastq/{file} | "
        f"samtools sort -o {file}.bam && "
        f"samtools index {file}.bam"
    )
    subprocess.run(hisat2_exec, shell=True)

    end_time = time.time()
    duration = end_time - start_time
    log_entry = f"{timestamp} - Finished processing {sample_name} in {duration:.2f} seconds\n"
    with open(logfile, "a") as log:
        log.write(log_entry)
    print(f"Processed {sample_name} in {duration:.2f} seconds")


print("All commands executed successfully")
```

Quantifying reads using FeatureCounts
-------------------------------------
The reads are quantified using FeatureCounts, which uses genomic coordinates in BAM files and maps to genomic features based on annotation in the reference genome.
```bash
featureCounts -s 0 -a ../fastq/Homo_sapiens.GRCh38.114.gtf \
        -o ../fastq/{bam}_featurecounts.txt \
        {bam}"
```
` -s ` : Indicate if strand-specific read counting should be performed.
` -a ` : The annotation file to be used
` -o ` : Name of the output directory