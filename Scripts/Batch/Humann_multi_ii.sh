#!/bin/bash

# Define input, output, and job directories
input_dir="/n/holylfs05/LABS/nguyen_lab/Everyone/users/zhijih/Madagascar/biobakery_output_new/rainforest/kneaddata/main/"
output_dir="/n/holystore01/LABS/huttenhower_lab/Users/Dongyu/PET_Projetc/Results/Madagascar_human/"
job_dir="./humann_jobs"
sample_list="/n/netscratch/huttenhower_lab/Lab/Users/Dongyu/Madagascar/List.csv"  # CSV file containing sample names (one per line)

# Create job directory if it doesn't exist
mkdir -p "$job_dir"

# Check if sample list file exists
if [[ ! -f "$sample_list" ]]; then
    echo "Error: Sample list file '$sample_list' not found!"
    exit 1
fi

# Read sample names from the list and process matching files
while IFS= read -r sample_name; do

    # Construct the expected file path based on the sample name
    file="${input_dir}${sample_name}.fastq.gz"

    # Check if the FASTQ file exists for the sample
    if [[ ! -f "$file" ]]; then
        echo "Warning: File $file not found, skipping..."
        continue
    fi

    # Define job script name
    job_script="${job_dir}/${sample_name}.slurm"

    # Ensure output directory for sample exists
    mkdir -p "${output_dir}${sample_name}"

    # Create SLURM job script for each valid file
    cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-12:00
#SBATCH -p huttenhower
#SBATCH --mem=90G
#SBATCH -o "${output_dir}${sample_name}/${sample_name}_%j.out"
#SBATCH -e "${output_dir}${sample_name}/${sample_name}_%j.err"
#SBATCH --mail-user=dongyu_wang@hsph.harvard.edu
#SBATCH --mail-type=ALL

# Load environment and required modules
source /n/huttenhower_lab/tools/hutlab/src/hutlabrc_rocky8.sh
hutlab load rocky8/biobakery_workflows/3.1.0-devel-dependsUpdate
hutlab load rocky8/metaphlan4/4.0.6_vOct22_fixed
hutlab load rocky8/humann4/4.0-alpha-1-final

# Check if HUMAnN is available before running
if ! command -v humann &> /dev/null; then
    echo "HUMAnN not found in the environment"
    exit 1
fi

# Run HUMAnN with the specified input and output
humann --threads 20 --input "$file" --output "${output_dir}${sample_name}/"

echo "HUMAnN analysis completed for $sample_name"
EOF

    # Submit the job to SLURM scheduler
    sbatch "$job_script"
    echo "Submitted job for: $sample_name"

done < "$sample_list"

echo "All jobs have been submitted successfully."
