import os
import zipfile
import re

def extract_fastqc_data(zip_path):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # Get a list of FastQC data files
        fastqc_files = [name for name in zip_ref.namelist() if re.match(r'.*/fastqc_data\.txt', name)]

        # Process each FastQC data file
        results = []
        for fastqc_file in fastqc_files:
            with zip_ref.open(fastqc_file) as file:
                # Read each FastQC data file
                fastqc_data = file.read().decode('utf-8')

                # Extract specific information from each FastQC data file
                total_sequences = re.search(r'Total Sequences\s+(\d+)', fastqc_data).group(1)
                gc_content = re.search(r'%GC\s+(\d+)', fastqc_data).group(1)
                deduplicated_percentage = re.search(r'Total Deduplicated Percentage\s+(\d+\.\d+)', fastqc_data).group(1)
                adapter_content = re.search(r'Adapter Content\s+(\S+)', fastqc_data).group(1)

                # Handle the case where 'Kmer Content' might not be present
                kmer_match = re.search(r'Kmer Content\s+(\S+)', fastqc_data)
                kmer_content = kmer_match.group(1) if kmer_match else "N/A"

                # Store the results
                results.append({
                    'Total Sequences': total_sequences,
                    '%GC': gc_content,
                    'Total Deduplicated Percentage': deduplicated_percentage,
                    'Adapter Content': adapter_content,
                    'Kmer Content': kmer_content
                })

        return results

def process_zip_folders(input_path, output_path):
    # Create the output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)

    # List all zip files in the input path
    zip_files = [f for f in os.listdir(input_path) if f.endswith('.zip')]

    # Process each zip file and extract FastQC data
    qc_summary = []
    for zip_file in zip_files:
        zip_path = os.path.join(input_path, zip_file)
        results = extract_fastqc_data(zip_path)

        # Write the summary to the output file
        qc_summary.append(f"Zip File: {zip_file}\n")
        for result in results:
            qc_summary.append(f"Total Sequences: {result['Total Sequences']}\n")
            qc_summary.append(f"%GC: {result['%GC']}\n")
            qc_summary.append(f"Total Deduplicated Percentage: {result['Total Deduplicated Percentage']}\n")
            qc_summary.append(f"Adapter Content: {result['Adapter Content']}\n")
            qc_summary.append(f"Kmer Content: {result['Kmer Content']}\n")
            qc_summary.append('-' * 50 + '\n')

    # Write the summary to qc_summary.txt
    summary_path = os.path.join(output_path, 'qc_summary.tsv')
    with open(summary_path, 'w') as summary_file:
        summary_file.writelines(qc_summary)

    print(f"Summary written to: {summary_path}")

if __name__ == "__main__":
    input_path = "home/anum/cell_lines_analysisoutput"
    output_path = "home/anum/cell_lines_analysisoutput/FastQC"

    process_zip_folders(input_path, output_path)