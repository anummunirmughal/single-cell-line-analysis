# single cell line analysis
 This repository contains codes for variant calling and filtering for single cell line exome and genome sequences. variant_calling.sh script can analyse single fastq file only.

 Step 1:perform fastqc analysis
 Step 2: Run parse_FastQC.py to determine whether your fastq files need trimmomatic analysis or not

After running the command analyse the output of parse_FastQC 
#Adapters: If FastQC reports fail for adapter content, trimming is necessary.
#Low Quality Bases: If FastQC reports low quality in the ends of reads, you may need to trim low-quality bases using tools like Trimmomatic.
#Overrepresented Sequences: If FastQC reports overrepresented sequences (other than adapters), trimming may be needed.
#GC content: if FastQC reports low GC content, trimming is necessary.
