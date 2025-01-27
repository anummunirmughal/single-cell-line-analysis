#!/bin/bash
#SBATCH --partition=biosoc2
#SBATCH --time=UNLIMITED
#SBATCH --cpus-per-task=20
#SBATCH --array=1
#SBATCH --output=1output.%A.%a.out
#SBATCH --error=1error.%A.%a.err

### Single-end Sequencing Analysis Script ###

# Paths
SOFTWARE_DIR="/home/anum/software"
INPUT_DIR="/home/anum/cell_lines_analysis/TNBC"
OUTPUT_DIR="/home/anum/cell_lines_analysisoutput"

DB_PATH="$SOFTWARE_DIR/datasets/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf"
FASTQC_PATH="$SOFTWARE_DIR/fastqc/FastQC"
BWA_PATH="$SOFTWARE_DIR/bwa/bwa-0.7.17"
SAMTOOLS_PATH="$SOFTWARE_DIR/samtools/samtools-1.15.1"
PICARD_PATH="$SOFTWARE_DIR/picard_tools/picard"
GATK_PATH="$SOFTWARE_DIR/gatk/gatk-4.2.6.1/gatk"
SNPEFF_PATH="$SOFTWARE_DIR/snpeff/snpEff"
GENOME="$SOFTWARE_DIR/datasets/hg38/Homo_sapiens_assembly38.fasta"

# Parameters
CPU_NUM=4
VEP_CPU_NUM=4

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Main workflow
echo "Starting the workflow..."

# Step 1: FastQC
echo "Step 1: Running FastQC..."
"$FASTQC_PATH/fastqc" -o "$OUTPUT_DIR" "$INPUT_DIR/SRR28563511.fastq"
if [ $? -ne 0 ]; then
    echo "Error: FastQC failed!"
    exit 1
fi

# Step 2: Align reads with BWA
echo "Step 2: Aligning reads with BWA..."
"$BWA_PATH/bwa" mem -t "$CPU_NUM" "$GENOME" "$INPUT_DIR/SRR28563511.fastq" > "$OUTPUT_DIR/aligned_reads.sam"
if [ $? -ne 0 ]; then
    echo "Error: BWA alignment failed!"
    exit 1
fi

# Step 3: Sort SAM to BAM
echo "Step 3: Sorting BAM file..."
"$SAMTOOLS_PATH/samtools" sort -@ "$CPU_NUM" -o "$OUTPUT_DIR/sorted_reads.bam" "$OUTPUT_DIR/aligned_reads.sam"
if [ $? -ne 0 ]; then
    echo "Error: BAM sorting failed!"
    exit 1
fi

# Step 4: Index BAM file
echo "Step 4: Indexing BAM file..."
"$SAMTOOLS_PATH/samtools" index "$OUTPUT_DIR/sorted_reads.bam"
if [ $? -ne 0 ]; then
    echo "Error: BAM indexing failed!"
    exit 1
fi

# Step 5: FixMateInformation
echo "Running Picard FixMateInformation..."
java -jar "$PICARD_PATH/picard.jar" FixMateInformation \
  -I "$OUTPUT_DIR/sorted_reads.bam" \
  -O "$OUTPUT_DIR/aligned_reads_fixmate.bam" \
  --TMP_DIR "$OUTPUT_DIR/temp" \
  --SORT_ORDER coordinate \
  --CREATE_INDEX true

# Step 6: Validate BAM file
echo "Validating BAM file..."
java -jar "$PICARD_PATH/picard.jar" ValidateSamFile \
  -I "$OUTPUT_DIR/aligned_reads_fixmate.bam" \
  -O "$OUTPUT_DIR/aligned_reads_validate_report.txt" \
  --MODE SUMMARY \
  --INDEX_VALIDATION_STRINGENCY EXHAUSTIVE

# Step 7: MarkDuplicates
echo "Marking duplicates..."
java -jar "$PICARD_PATH/picard.jar" MarkDuplicates \
  -I "$OUTPUT_DIR/aligned_reads_fixmate.bam" \
  -O "$OUTPUT_DIR/aligned_reads_dedup.bam" \
  --METRICS_FILE "$OUTPUT_DIR//dedup_metrics.txt" \
  --REMOVE_DUPLICATES false \
  --CREATE_INDEX true \
  --TMP_DIR "$OUTPUT_DIR/temp"

# Step 8: Index final BAM file
echo "Indexing the final BAM file..."
"$SAMTOOLS_PATH/samtools" index -@ "$THREADS" "$OUTPUT_DIR/aligned_reads_dedup.bam"