#!/bin/bash

# Define input, output, and job directories
input_dir="/n/holystore01/LABS/huttenhower_lab/Lab/data/hmp1_II_Stool_concat_postknead/"
output_dir="/n/netscratch/huttenhower_lab/Lab/Users/Dongyu/HMP1-II_humann/"
job_dir="/n/netscratch/huttenhower_lab/Lab/Users/Dongyu/humann_jobs"
file_list="/n/netscratch/huttenhower_lab/Lab/Users/Dongyu/HMP1-II/List.csv"  # List of specific sample names

# Ensure job directory is created and accessible
if ! mkdir -p "$job_dir"; then
    echo "Error: Cannot create directory $job_dir. Check permissions!" >&2
    exit 1
fi

batch_size=5  # Number of jobs per batch
count=0       # Counter for tracking number of samples
batch=()      # Array to store batch samples

# Read file names line by line from the list
while IFS= read -r sample || [[ -n "$sample" ]]; do
    file="${input_dir}${sample}.fastq.gz"

    # Check if the FASTQ file exists before proceeding
    if [[ -f "$file" ]]; then
        batch+=("$sample")
        ((count++))
    else
        echo "Warning: File $file not found. Skipping..."
        continue
    fi

    # If batch size is reached, submit jobs and wait
    if (( count % batch_size == 0 )); then
        echo "Submitting batch of $batch_size jobs..."
        job_ids=()  # Store submitted job IDs

        for sample in "${batch[@]}"; do
            job_script="${job_dir}/${sample}.slurm"

            # Ensure output directory exists
            mkdir -p "${output_dir}${sample}/"

            # Create SLURM job script
            cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-12:00
#SBATCH -p test
#SBATCH --mem=100G
#SBATCH -o ${output_dir}${sample}/Humann_${sample}.out
#SBATCH -e ${output_dir}${sample}/Humann_${sample}.err
#SBATCH --mail-user=dongyu_wang@hsph.harvard.edu
#SBATCH --mail-type=ALL

# Load environment and required modules
source /n/huttenhower_lab/tools/hutlab/src/hutlabrc_rocky8.sh
hutlab load rocky8/biobakery_workflows/3.1.0-devel-dependsUpdate
hutlab load rocky8/metaphlan4/4.0.6_vOct22_fixed
hutlab load rocky8/humann4/4.0-alpha-1-final

# Run HUMAnN with the specified input and output
humann --threads 20 --input "$file" --output "${output_dir}${sample}/"

echo "HUMAnN analysis completed for $sample"
EOF

            # Submit job and capture job ID
            job_id=$(sbatch "$job_script" | awk '{print $4}')
            job_ids+=("$job_id")
            echo "Submitted job for: $sample (Job ID: $job_id)"
        done

        # Wait for the batch to finish before proceeding
        echo "Waiting for batch to complete..."
        for jid in "${job_ids[@]}"; do
            while squeue -j "$jid" &> /dev/null; do
                sleep 60  # Check every 60 seconds
            done
        done

        echo "Batch completed. Proceeding to the next batch..."

        # Reset batch array
        batch=()
    fi
done < "$file_list"

# Submit remaining jobs if any
if [ ${#batch[@]} -gt 0 ]; then
    echo "Submitting final batch of ${#batch[@]} jobs..."
    job_ids=()

    for sample in "${batch[@]}"; do
        job_script="${job_dir}/${sample}.slurm"
        mkdir -p "${output_dir}${sample}/"

        # Create SLURM job script
        cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-12:00
#SBATCH -p test
#SBATCH --mem=100G
#SBATCH -o ${output_dir}${sample}/Humann_${sample}.out
#SBATCH -e ${output_dir}${sample}/Humann_${sample}.err
#SBATCH --mail-user=dongyu_wang@hsph.harvard.edu
#SBATCH --mail-type=ALL

# Load environment and required modules
source /n/huttenhower_lab/tools/hutlab/src/hutlabrc_rocky8.sh
hutlab load rocky8/biobakery_workflows/3.1.0-devel-dependsUpdate
hutlab load rocky8/metaphlan4/4.0.6_vOct22_fixed
hutlab load rocky8/humann4/4.0-alpha-1-final

# Run HUMAnN with the specified input and output
humann --threads 20 --input "$file" --output "${output_dir}${sample}/"

echo "HUMAnN analysis completed for $sample"
EOF

        # Submit job and capture job ID
        job_id=$(sbatch "$job_script" | awk '{print $4}')
        job_ids+=("$job_id")
        echo "Submitted job for: $sample (Job ID: $job_id)"
    done

    echo "Waiting for final batch to complete..."
    for jid in "${job_ids[@]}"; do
        while squeue -j "$jid" &> /dev/null; do
            sleep 60  # Check every 60 seconds
        done
    done

    echo "Final batch completed."
fi

echo "All jobs have been submitted successfully."
