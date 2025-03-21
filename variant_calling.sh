#!/bin/bash
<<<<<<< Updated upstream
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
BCFTOOLS_PATH="$SOFTWARE_DIR/bcftools-1.19/bin/bcftools"

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
"$SAMTOOLS_PATH/samtools" index -@ "$CPU_NUM" "$OUTPUT_DIR/aligned_reads_dedup.bam"

#Step 9: Basecalibration
echo "Step 9: basecalibration..."
"$GATK_PATH" BaseRecalibrator \
  -I "$OUTPUT_DIR/aligned_reads_dedup.bam" \
  -R "$GENOME" \
  --known-sites /shared/teams/hpc-wasslab/updated_pipeline/data/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites /shared/teams/hpc-wasslab/updated_pipeline/data/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites /shared/teams/hpc-wasslab/updated_pipeline/data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites /shared/teams/hpc-wasslab/updated_pipeline/data/common_all_20180418.vcf \
  -O "$OUTPUT_DIR/recal_chr1.table" 

echo "Completed Step 1 Base Recalibrator"

# Step 10: Apply the model to adjust the base quality scores
echo "Step 10: Scoring of calibration"
"$GATK_PATH" ApplyBQSR \
  -I "$OUTPUT_DIR/aligned_reads_dedup.bam" \
  -R "$GENOME" \
  --bqsr-recal-file "$OUTPUT_DIR/recal_chr1.table" \
  -O "$OUTPUT_DIR/BQSR/chr1_sort_dup_bqsr.bam" \
if [ $? -ne 0 ]; then
    echo "Error: ApplyBQSR failed!"
    exit 1
fi  
  
# Step 11: Generate VCF using bcftools mpileup and call variants
echo "Step 11: Running bcftools mpileup and variant calling..."
"$BCFTOOLS_PATH" mpileup -Ou -f "$GENOME" "$OUTPUT_DIR/chr1_sort_dup_bqsr.bam" | \
  "$BCFTOOLS_PATH" call -mv -Ob -o "$OUTPUT_DIR/raw_calls.bcf"
if [ $? -ne 0 ]; then
    echo "Error: Variant calling failed!"
    exit 1
fi

# Step 12: Filter variants based on quality
echo "Step 12: Filtering variants..."
"$BCFTOOLS_PATH" filter -i '%QUAL>=30' \
  -Ob -o "$OUTPUT_DIR/filtered_calls.bcf" "$OUTPUT_DIR/raw_calls.bcf"

# Step 13: Convert to VCF format
echo "Step 13: Converting to VCF format..."
"$BCFTOOLS_PATH" view "$OUTPUT_DIR/filtered_calls.bcf" > "$OUTPUT_DIR/final_variants.vcf"
if [ $? -ne 0 ]; then
    echo "Error: Filtering or VCF conversion failed!"
    exit 1
fi

echo "Variant calling and filtering completed!"
  
=======

#SBATCH --time=UNLIMITED
#SBATCH --cpus-per-task=20
#SBATCH --array=1-4 ### Change this for the total samples to analyze
#SBATCH --output=*.%A.%a.out

### sequencing Project - Script 01 ###

# Paths
SOFTWARE_DIR="/home/anum/software"
INPUT_DIR="/home/anum/single_cell_lines/input"
OUTPUT_DIR="/home/anum/single_cell_lines/output"

# Default paths for software tools
INDEX_PATH="$SOFTWARE_DIR/hg38/Homo_sapiens_assembly38"
DB_PATH="$SOFTWARE_DIR/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf"
FASTQC_PATH="$SOFTWARE_DIR/FastQC"
TRIMMOMATIC_PATH="$SOFTWARE_DIR/Trimmomatic-0.39"
BWA_PATH="$SOFTWARE_DIR/bwa-0.7.17"
PICARD_PATH="$SOFTWARE_DIR/picard-2.17.10"
GATK_PATH="$SOFTWARE_DIR/gatk-4.2.6.1"
SAMTOOLS_PATH="$SOFTWARE_DIR/samtools-1.15.1"
BCFTOOLS_PATH="$SOFTWARE_DIR/bcftools-1.15.1"
VEP_PATH="$SOFTWARE_DIR/ensembl-vep"

# Parameters
REMOVE_INTERMEDIATE="n"
REMOVE_LOGS="n"
PRINT_COMMANDS="n"

# Usage information
USAGE="
Usage: $(basename "$0") -i INPUT_DIR -o OUTPUT_DIR [options]

Required Arguments:
  -i  Input directory containing FASTQ files (default: $INPUT_DIR)
  -o  Output directory for results (default: $OUTPUT_DIR)

Optional Arguments:
  -rmf   Remove intermediate files (y/n) [default: $REMOVE_INTERMEDIATE]
  -rme   Remove error logs (y/n) [default: $REMOVE_LOGS]
  -cpu   Number of CPUs to use in parallel (default: $CPU_NUM)
  -vcpu  Number of CPUs for VEP parallelization (default: $VEP_CPU_NUM)
  -cl    Print commands only (y/n) [default: $PRINT_COMMANDS]
"


# Print settings
echo "Input Directory: $INPUT_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo "Remove Intermediate Files: $REMOVE_INTERMEDIATE"
echo "Remove Logs: $REMOVE_LOGS"
echo "CPUs: $CPU_NUM"
echo "VEP CPUs: $VEP_CPU_NUM"
echo "Commands Only: $PRINT_COMMANDS"

# Main workflow
if [[ "$PRINT_COMMANDS" == "y" ]]; then
    echo "Commands will be printed but not executed."
else
    echo "Running the workflow..."
    
    # Step 1: FastQC
    echo "Step 1: Running FastQC..."
    "$FASTQC_PATH/fastqc" -o "$OUTPUT_DIR" "$INPUT_DIR"/*_R1.fastq "$INPUT_DIR"/*_R2.fastq
    
    # Step 2: Trimmomatic
    echo "Step 2: Running Trimmomatic..."
    java -jar "$TRIMMOMATIC_PATH/trimmomatic.jar" PE -threads "$CPU_NUM" \
        "$INPUT_DIR"/*_R1.fastq "$INPUT_DIR"/*_R2.fastq \
        "$OUTPUT_DIR"/trimmed_R1.fastq "$OUTPUT_DIR"/unpaired_R1.fastq \
        "$OUTPUT_DIR"/trimmed_R2.fastq "$OUTPUT_DIR"/unpaired_R2.fastq \
        ILLUMINACLIP:"$TRIMMOMATIC_PATH/adapters/TruSeq3-PE.fa":2:30:10

    # Step 3: Align reads with BWA
    echo "Step 3: Aligning reads with BWA..."
    "$BWA_PATH/bwa" mem -t "$CPU_NUM" "$INDEX_PATH" \
        "$OUTPUT_DIR"/trimmed_R1.fastq "$OUTPUT_DIR"/trimmed_R2.fastq > "$OUTPUT_DIR"/aligned_reads.sam

    # Step 4: Convert SAM to BAM
    echo "Step 4: Converting SAM to BAM..."
    "$SAMTOOLS_PATH/samtools" view -bS "$OUTPUT_DIR"/aligned_reads.sam > "$OUTPUT_DIR"/aligned_reads.bam

    # Step 5: Sort BAM file
    echo "Step 5: Sorting BAM file..."
    "$SAMTOOLS_PATH/samtools" sort -@ "$CPU_NUM" -o "$OUTPUT_DIR"/sorted_reads.bam "$OUTPUT_DIR"/aligned_reads.bam

    # Step 6: Index BAM file
    echo "Step 6: Indexing BAM file..."
    "$SAMTOOLS_PATH/samtools" index "$OUTPUT_DIR"/sorted_reads.bam

    # Step 7: Mark duplicates with Picard
    echo "Step 7: Marking duplicates with Picard..."
    java -jar "$PICARD_PATH/picard.jar" MarkDuplicates \
        I="$OUTPUT_DIR"/sorted_reads.bam \
        O="$OUTPUT_DIR"/dedup_reads.bam \
        M="$OUTPUT_DIR"/metrics.txt

    # Step 8: Index deduplicated BAM
    echo "Step 8: Indexing deduplicated BAM file..."
    "$SAMTOOLS_PATH/samtools" index "$OUTPUT_DIR"/dedup_reads.bam

    # Step 9: Print BAM statistics
    echo "Step 9: Generating BAM statistics..."
    "$SAMTOOLS_PATH/samtools" flagstat "$OUTPUT_DIR"/dedup_reads.bam > "$OUTPUT_DIR"/bam_statistics.txt

    # Step 10: Coverage analysis
    echo "Step 10: Calculating coverage..."
    "$SAMTOOLS_PATH/samtools" depth "$OUTPUT_DIR"/dedup_reads.bam > "$OUTPUT_DIR"/coverage.txt

    # Step 11: Variant calling with GATK
    echo "Step 11: Variant calling with GATK..."
    java -jar "$GATK_PATH/gatk.jar" HaplotypeCaller \
        -R "$INDEX_PATH" \
        -I "$OUTPUT_DIR"/dedup_reads.bam \
        -O "$OUTPUT_DIR"/variants.vcf

    # Step 12: Annotate variants with VEP
    echo "Step 12: Annotating variants with VEP..."
    perl "$VEP_PATH/vep" -i "$OUTPUT_DIR"/variants.vcf -o "$OUTPUT_DIR"/annotated_variants.vcf \
        --cache --dir_cache "$VEP_PATH/cache" --fork "$VEP_CPU_NUM"
fi

# Cleanup
if [[ "$REMOVE_INTERMEDIATE" == "y" ]]; then
    echo "Removing intermediate files..."
    rm -f "$OUTPUT_DIR"/*.sam "$OUTPUT_DIR"/*.bam "$OUTPUT_DIR"/*.txt
fi

if [[ "$REMOVE_LOGS" == "y" ]]; then
    echo "Removing log files..."
    rm -f "$OUTPUT_DIR"/*.log
fi

echo "Workflow completed."
>>>>>>> Stashed changes
