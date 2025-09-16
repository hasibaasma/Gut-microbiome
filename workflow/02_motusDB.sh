#!/bin/bash
# ============================================================
# 02_motusdb.sh
# Taxonomic profiling of pre-trimmed FASTQ files using mOTUs
# ============================================================

set -euo pipefail

# -----------------------------
# User settings
# -----------------------------
THREADS=8
TRIMMED_DIR="samples_trimmed"   # directory with *_trimmed.fastq.gz
PROFILE_DIR="${TRIMMED_DIR}/profiles"
MERGED_TABLE="merged_mOTUs_8_table.tsv"

# -----------------------------
# 1. Install mOTUs (only once)
#    NOTE: This requires Python ≥3.8 and pip
# -----------------------------
echo ">>> Installing mOTUs via pip (skip if already installed) ..."
pip install --upgrade motu-profiler

# -----------------------------
# 2. Download mOTUs database (only once)
# -----------------------------
echo ">>> Downloading mOTUs database ..."
motus downloadDB

# -----------------------------
# 3. Profile each sample
# -----------------------------
echo ">>> Running mOTUs profiling on all trimmed FASTQ files ..."
mkdir -p "${PROFILE_DIR}"

for fq in "${TRIMMED_DIR}"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    echo "---- Profiling sample: ${sample}"
    motus profile -s "$fq" -n "$sample" -o "${PROFILE_DIR}/${sample}.motus.tsv" -t "${THREADS}"
done

# -----------------------------
# 4. Merge all per-sample tables
# -----------------------------
echo ">>> Merging mOTUs profiles into a single abundance table ..."
cd "${PROFILE_DIR}"
motus merge -i *.motus.tsv -o "../${MERGED_TABLE}"
cd - > /dev/null

echo ">>> Done."
echo "Merged table available at: ${TRIMMED_DIR}/${MERGED_TABLE}"
#!/bin/bash
# ============================================================
# 02_motusdb.sh
# Taxonomic profiling of pre-trimmed FASTQ files using mOTUs
# ============================================================

set -euo pipefail

# -----------------------------
# User settings
# -----------------------------
THREADS=8
TRIMMED_DIR="samples_trimmed"   # directory with *_trimmed.fastq.gz
PROFILE_DIR="${TRIMMED_DIR}/profiles"
MERGED_TABLE="merged_mOTUs_8_table.tsv"

# -----------------------------
# 1. Install mOTUs (only once)
#    NOTE: This requires Python ≥3.8 and pip
# -----------------------------
echo ">>> Installing mOTUs via pip (skip if already installed) ..."
pip install --upgrade motu-profiler

# -----------------------------
# 2. Download mOTUs database (only once)
# -----------------------------
echo ">>> Downloading mOTUs database ..."
motus downloadDB

# -----------------------------
# 3. Profile each sample
# -----------------------------
echo ">>> Running mOTUs profiling on all trimmed FASTQ files ..."
mkdir -p "${PROFILE_DIR}"

for fq in "${TRIMMED_DIR}"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    echo "---- Profiling sample: ${sample}"
    motus profile -s "$fq" -n "$sample" -o "${PROFILE_DIR}/${sample}.motus.tsv" -t "${THREADS}"
done

# -----------------------------
# 4. Merge all per-sample tables
# -----------------------------
echo ">>> Merging mOTUs profiles into a single abundance table ..."
cd "${PROFILE_DIR}"
motus merge -i *.motus.tsv -o "../${MERGED_TABLE}"
cd - > /dev/null

echo ">>> Done."
echo "Merged table available at: ${TRIMMED_DIR}/${MERGED_TABLE}"
