#Go to your aligned reads folder
import os
import subprocess
import time
from datetime import datetime 

# dwn = f"wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz"
# unzip = f"gunzip Homo_sapiens.GRCh38.114.gtf.gz"
# subprocess.run(dwn, shell=True)
# subprocess.run(unzip, shell=True)

bam_files = [
    # "LNCAP_Normoxia_S1.bam",
    # "LNCAP_Normoxia_S2.bam",
    # "LNCAP_Hypoxia_S1.bam",
    # "LNCAP_Hypoxia_S2.bam",
    # "PC3_Normoxia_S1.bam",
    "PC3_Normoxia_S2.bam",
    "PC3_Hypoxia_S1.bam",
    "PC3_Hypoxia_S2.bam"
]

# Loop over all BAM files
for bam in bam_files:
    start=time.time()  # start time

    print(f"{start} Processing {bam} ...")
    featurecounts = f"featureCounts -s 0 -a /mnt/d/SMau/bulkrnaseq_tutorial/fastq/Homo_sapiens.GRCh38.114.gtf \
        -o /mnt/d/SMau/bulkrnaseq_tutorial/fastq/{bam}_featurecounts.txt \
        {bam}"
    subprocess.run(featurecounts, shell=True)

    end=time.time()  # end time
    runtime=(( (end - start) / 60 ))  # in minutes

    print(f"✅ Completed {bam} in {runtime} minutes.")
    