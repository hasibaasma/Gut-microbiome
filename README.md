

#  README

## 1. Data sources

This task is based on publicly available data from a study of **Dynamics of Gut Microbiota After Fecal Microbiota Transplantation in Ulcerative Colitis: Success Linked to Control of Prevotellaceae**. (https://pmc.ncbi.nlm.nih.gov/articles/PMC11836888/#s27) The dataset includes multiple samples under different conditions (209 samples) and was originally sequenced using **Illumina NovaSeq platform (100 bp single end read)**.
The subsampled and cleaned FASTQs are stored in `data/` and are used as the inputs for the workflow.

---

## 2. How to download

The raw data can be downloaded from NCBI’s SRA database using the SRA Toolkit. For example, to download a sample (SRR27827162) into the data/ folder:
```
mkdir -p data

# Download FASTQ using fasterq-dump
fasterq-dump SRR27827162 \
  --threads 8 \
  --outdir data

```
This will generate
```
data/SRR27827162.fastq
```

---

## 3. Pre-processing / subsampling

To make the dataset smaller and faster to process for testing, I randomly subsampled reads using seqtk. So to keep 10% of reads this was done  
```.
for f in data/*.fastq; do
    base=$(basename $f .fastq)
    seqtk sample -s100 $f 0.1 > data/${base}_sub.fastq
done

```

---

## 4. How the workflow works
Removal of host genome uses **Bowtie2** to map reads to reference genome  
Filter unmapped paired reads from output using **Samtools**  
Quality Trimming, adapter removal and low complexity filtering is done with **\ Fastp** 
Reports of mapping and quality trimming are generated with **\ MultiQC**  
The input files is stored in `data/`.  

---
### Step 1 - Organize your files
**Purpose:** Ensure raw single-end sequencing reads are consistently named for downstream workflow compatibility.  
**Tools:** None (manual or scripting via mv, rename, or bash loop).  
**Inputs:** Raw FASTQ files with arbitrary names from sequencing provider.  
**Outputs:** Renamed FASTQ files following convention:
sample1.fastq.gz
sample2.fastq.gz
**Command:**
```
# Example using rename to standardize to sampleX.fastq.gz
rename 's/(.*)\.fastq\.gz$/sample_$1.fastq.gz/' *.fastq.gz 
```
---

### Step 2 - Build host genome index

**Purpose:** Create a searchable index of the host reference genome for efficient read alignment (needed for host read removal).  
**Tools:** Bowtie2   
**Inputs:** Host reference genome FASTA file (e.g., GRCh38.fa).   
**Outputs:** Bowtie2 index files (host_reference.1.bt2, host_reference.2.bt2, …, host_reference.rev.2.bt2).  
 
**Command:**

```
bowtie2-build GRCh38.fa host_reference
```

---
### Step 3 – Quality Trimming & Adapter Removal
 
**Purpose:** Improve read quality by removing low-quality bases, adapters, and very short reads. This ensures better mapping and reduces false positives in downstream analyses.  
**Tools:** fastp  
**Inputs:**  
Paired-end FASTQ files after host removal  
  sample1_hostRemoved_R1.fastq.gz
  sample1_hostRemoved_R2.fastq.gz  
**Outputs:**  
Cleaned FASTQ files  
  sample1_trimmed_R1.fastq.gz  
  sample1_trimmed_R2.fastq.gz  
QC reports    
  sample1_fastp.html (interactive QC report)
  sample1_fastp.json (machine-readable QC summary)  
**Command:**
```
fastp \
  --in1 sample_hostRemoved_R1.fastq.gz \
  --in2 sample_hostRemoved_R2.fastq.gz \
  --out1 sample_trimmed_R1.fastq.gz \
  --out2 sample_trimmed_R2.fastq.gz \
  --cut_right --cut_window_size 4 --cut_mean_quality 20 \
  -l 50 \
  --detect_adapter_for_pe \
  -y \
  --thread 8 \
  --html sample_fastp.html \
  --json sample_fastp.json

```

---

### Step 4 - Generate Summary Report

**Purpose:** Aggregate QC and alignment results from multiple tools (Bowtie2, Samtools, Fastp) into a single, interactive report for easier interpretation and comparison across samples.  
**Tools:** MultiQC  
**Inputs:**  
Log files and reports from previous steps, e.g.:
  Bowtie2 alignment logs (*.log)
  Samtools filtering logs (if any)
  Fastp QC outputs (*.html, *.json)
**Outputs:**  
  Combined QC report:  
  multiqc_report.html (interactive summary)  
  multiqc_data/ (supporting data files)  
**Command:**  
```
multiqc . -o multiqc_report
```

All of the above steps can be executed in one go using following command. Make sure you execute this in data folder.
```
bash workflow.sh
```
---


