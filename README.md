# Gut Microbiome mOTUs Workflow
## 1️. Study & Data Sources
This workflow reproduces analyses from the study  
“Dynamics of Gut Microbiota After Fecal Microbiota Transplantation in Ulcerative Colitis: Success Linked to Control of Prevotellaceae”
(PMC11836888).  
•	Sequencing platform: Illumina NovaSeq (100 bp single-end)  
•	Dataset: 209 publicly available gut‐microbiome samples.  
•	Inputs: Subsampled and cleaned FASTQ files in data/.  
 
## 2️. Download Raw Data
Use the SRA Toolkit to fetch FASTQ files from NCBI SRA:  
```
mkdir -p data
fasterq-dump SRR27827162 --threads 8 --outdir data
```
This produces:
```
data/SRR27827162.fastq
 ```
## 3. Subsampling
To reduce runtime, you may down-sample reads to 10 % with seqtk:  
```
for f in data/*.fastq; do
    base=$(basename "$f" .fastq)
    seqtk sample -s100 "$f" 0.1 > "data/${base}_sub.fastq"
done
``` 
## 4️. Workflow Overview
The complete analysis consists of three main stages:  
1.	Pre-processing: Host-read removal, quality control and summary reports.  
2.	Taxonomic profiling: mOTUs v3 for species-level abundance tables.  
3.	Downstream R analysis: Microbiota composition, PCA ordination and diversity metrics.  
All commands below assume a Linux/macOS environment with required tools installed.  
 
### 4.1 Pre-processing
#### a. Organize FASTQ Files  
Rename or reorganize raw reads into a consistent pattern such as:  
```
sample1.fastq.gz
sample2.fastq.gz
```
(to simplify automation).     
#### b. Build Host Genome Index   
Create a Bowtie2 index for the host genome (e.g. human GRCh38):    
```
bowtie2-build GRCh38.fa host_reference
```
#### c. Host Read Removal, Trimming & QC 
•	Bowtie2: map reads to host reference and keep unmapped reads.  
•	Samtools: extract unmapped reads.  
•	fastp: trim adapters, low-quality ends and short reads.  
•	MultiQC: generate a combined QC report.  
Example fastp command:  
```
fastp \
  --in1 sample_hostRemoved_R1.fastq.gz \
  --out1 sample_trimmed_R1.fastq.gz \
  --cut_right --cut_window_size 4 --cut_mean_quality 20 \
  -l 50 \
  --html sample_fastp.html \
  --json sample_fastp.json \
  --thread 8
```
Summarize results:
```
multiqc . -o multiqc_report
```
You can automate the entire pre-processing phase:
```
bash workflow.sh
```
Outputs: cleaned FASTQ files in samples_trimmed/, plus an interactive report multiqc_report.html.
 
### 4.2 mOTUs Taxonomic Profiling
#### Installation (first time only)
```
pip install --upgrade motu-profiler
motus downloadDB
```
#### Run profiling and merge tables
```
bash scripts/02_motusdb.sh
```
#### Outputs
• Individual profiles: samples_trimmed/profiles/*.motus.tsv
• Merged abundance table: samples_trimmed/merged_mOTUs_8_table.tsv
(rows = mOTUs, columns = samples)
This merged table is the input for downstream R analysis.
 
### 4.3 Post-processing and R-based Analysis
Convert the merged table into a phyloseq object and perform key ecological analyses:
•	Microbiota community composition – visualize relative abundances of major families.
•	Principal Component Analysis (PCA) – ordination of samples using Aitchison (CLR) distance to explore clustering patterns (e.g. donors vs. responders).
•	Simpson dominance – a diversity index reflecting whether a few taxa dominate or the community is evenly distributed.
#### Requirements
R ≥ 4.2 and its packages:

```
install.packages(c("phyloseq","ggplot2","vegan"))
if (!requireNamespace("microbiome", quietly = TRUE))
    BiocManager::install("microbiome")
```
#### Run analysis
```
Rscript scripts/03_postprocess.R
```
#### Outputs
• mOTUs_phyloseq.rds – phyloseq object for further custom analysis
• mOTUs_family_barplot.png – stacked barplot of relative abundances
• mOTUs_PCA.png – PCA ordination plot
• mOTUs_simpson_diversity.csv – Simpson dominance values per sample
These deliver a concise overview of the gut-microbiota community structure and diversity.
You can continue exploring differential abundance or longitudinal trends using the saved mOTUs_phyloseq.rds.

