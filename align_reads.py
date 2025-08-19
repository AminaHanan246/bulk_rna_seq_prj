import subprocess
import os
import time
from datetime import datetime

#Summarise report from multiple fastqc reports
multiqc_exec = f"multiqc fastqc_results/ -o multiqc_report/"
subprocess.run(multiqc_exec, shell=True)

#Running trimmomatic for one file and fastqc report
trim_exec = f"trimmomatic SE -threads 4 SRR7179504_pass.fastq.gz SRR7179504_trimmed.fastq.gz TRAILING:10 -phred33"
os.makedirs("fastqc_results", exist_ok=True)
fastqc_exec = f"fastqc SRR7179504_trimmed.fastq.gz -o fastqc_results/ --threads 4"
subprocess.run(trim_exec, shell=True)
subprocess.run(fastqc_exec, shell=True)

#Concatenating SRA files based on sample data
concat_exec = [
    f"cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz > LNCAP_Normoxia_S1.fastq.gz",
    f"cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz > LNCAP_Normoxia_S2.fastq.gz",
    f"cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz > LNCAP_Hypoxia_S1.fastq.gz",
    f"cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz > LNCAP_Hypoxia_S2.fastq.gz"
]

move_exec = [
    f"mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz",
    f"mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz",
    f"mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz",
    f"mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz"
]
for cmd in concat_exec:
    subprocess.run(cmd, shell=True)
for cmd in move_exec:
    subprocess.run(cmd, shell=True)

#downloading genome index file and unziping
download_exec = f"wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz"
unzip_exec = f"tar -xvzf grch38_genome.tar.gz"
subprocess.run(download_exec, shell=True)
subprocess.run(unzip_exec, shell=True)

files = [
    "LNCAP_Normoxia_S1.fastq.gz",
    "LNCAP_Normoxia_S2.fastq.gz",
    "LNCAP_Hypoxia_S1.fastq.gz",
    "LNCAP_Hypoxia_S2.fastq.gz",
    "PC3_Normoxia_S1.fastq.gz",
    "PC3_Normoxia_S2.fastq.gz",
    "PC3_Hypoxia_S1.fastq.gz",
    "PC3_Hypoxia_S2.fastq.gz"
]

#Alignment of reads via Hisat2
logfile = "alignment_log.txt"

for file in files:
    sample_name = os.path.basename(file).replace(".fastq.gz", "")
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"{timestamp} - Processing {sample_name}\n"
    print(log_entry)
    with open(logfile, "a") as log:
        log.write(log_entry)
    
    start_time = time.time()

    hisat2_exec = (
        f"hisat2 -p 6 -x grch38/genome -U {file} | "
        f"samtools sort -o {sample_name}.bam && "
        f"samtools index {sample_name}.bam"
    )
    print(f"Running command: hisat2")
    subprocess.run(hisat2_exec, shell=True)

    end_time = time.time()
    duration = end_time - start_time
    log_entry = f"{timestamp} - Finished processing {sample_name} in {duration:.2f} seconds\n"
    with open(logfile, "a") as log:
        log.write(log_entry)
    print(f"Processed {sample_name} in {duration:.2f} seconds")


print("All commands executed successfully")