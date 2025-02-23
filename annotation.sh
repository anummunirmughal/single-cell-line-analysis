#!/bin/bash
#SBATCH --partition=biosoc2
#SBATCH --time=UNLIMITED
#SBATCH --cpus-per-task=20
#SBATCH --array=1
#SBATCH --output=annotation_output.%A.%a.out
#SBATCH --error=1error.%A.%a.err

SOFTWARE_DIR="/home/anum/software"
WORKDIR="/home/anum/cell_lines_analysisoutput"
DB_PATH="${SOFTWARE_DIR}/datasets/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf"
data_PATH=""${SOFTWARE_DIR}/datasets"
SNPEFF_PATH="${SOFTWARE_DIR}/snpeff/snpEff"
GENOME="${SOFTWARE_DIR}/datasets/hg38/Homo_sapiens_assembly38.fasta"
BCFTOOLS_PATH="${SOFTWARE_DIR}/BCFTOOLS-1.19/bin/BCFTOOLS"
htslib_PATH="${SOFTWARE_DIR}/htslib-1.17"

mkdir -p "${WORKDIR}"/annotation/
java -Xmx4g -jar "${SNPEFF_PATH}" GRCh38.99 \
    -i vcf -o vcf \
    "${WORKDIR}"/Filter/final_variants.vcf \
    > "${WORKDIR}"/annotation/annotated_variants.vcf
echo "Annotation with snpEff completed"

#Fix: Compress & index the annotated file
"${htslib_PATH}"/bgzip -c "${WORKDIR}/annotation/annotated_variants.vcf" > "${WORKDIR}/annotation/annotated_variants.vcf.gz"
"${htslib_PATH}"/tabix -p vcf "${WORKDIR}/annotation/annotated_variants.vcf.gz"
echo "indexing finished"

# Annotate with gnomad
"${BCFTOOLS_PATH}" annotate -a "${data_PATH}"/gnomad/gnomad.genomes.v4.1.sites.chr1.vcf.bgz \
    -c INFO/AF "${WORKDIR}"/annotation/annotated_variants.vcf.gz | \

"${BCFTOOLS_PATH}" view -i 'INFO/AF<0.01' > "${WORKDIR}"/annotation/filtered_variants.vcf

# Sort the VCF file
"${BCFTOOLS_PATH}" sort -O v -o "${WORKDIR}/annotation/filtered_variants.sorted.vcf" "${WORKDIR}/annotation/filtered_variants.vcf"

# Compress and index the sorted VCF file
"${htslib_PATH}"/bgzip -c "${WORKDIR}/annotation/filtered_variants.sorted.vcf" > "${WORKDIR}/annotation/filtered_variants.vcf.gz"
"${htslib_PATH}"/tabix -p vcf "${WORKDIR}/annotation/filtered_variants.vcf.gz"

#Annotate with cosmic
"${BCFTOOLS_PATH}" annotate -a "${data_PATH}/cosmic/Cosmic_GenomeScreensMutant_v101_GRCh38.vcf.gz"  \
   -c ID,INFO "${WORKDIR}"/annotation/annotated_variants.vcf.gz > \
    "${WORKDIR}"/annotation/cosmic_annotated.vcf

# sort the VCF file
"${BCFTOOLS_PATH}" sort -O v -o "${WORKDIR}/annotation/cosmic_annotated.sorted.vcf" "${WORKDIR}/annotation/cosmic_annotated.vcf"

# Compress and index the sorted VCF file
"${htslib_PATH}"/bgzip -c "${WORKDIR}/annotation/cosmic_annotated.vcf" > "${WORKDIR}/annotation/cosmic_annotated.vcf.gz"
"${htslib_PATH}"/tabix -f vcf "${WORKDIR}/annotation/cosmic_annotated.vcf.gz"

# Concatenate both VCF files annotated with gnomad and cosmic
"${BCFTOOLS_PATH}" concat -a \
    "${WORKDIR}/annotation/filtered_variants.vcf.gz" \
    "${WORKDIR}/annotation/cosmic_annotated.vcf.gz" > "${WORKDIR}/annotation/final_filtered_variants.vcf"

# Compress & Index Final VCF
"${htslib_PATH}"/bgzip -c "${WORKDIR}/annotation/final_variants.vcf" > "${WORKDIR}/annotation/final_filtered_variants.vcf.gz"
"${htslib_PATH}"/tabix -p vcf "${WORKDIR}/annotation/final_filtered_variants.vcf.gz"

#Visualise the final variants file
"${BCFTOOLS_PATH}" view "${WORKDIR}"/annotation/final_filtered_variants.vcf.gz | less -S
