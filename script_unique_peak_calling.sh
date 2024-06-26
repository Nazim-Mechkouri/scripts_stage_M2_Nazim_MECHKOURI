#!/bin/bash
#SBATCH --job-name=peak_calling
#SBATCH --output=peak_calling_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

# Load modules
module load macs2
module load samtools

# Input parameters
BAM_FILE=$1
REFERENCE_GENOME_FILE=$2
OUTPUT_DIR=$3

# Create output directory if it does not exist
mkdir -p $OUTPUT_DIR

# Function to calculate genome size from reference genome
calculate_genome_size() {
    local genome_file=$1
    grep -v "^>" $genome_file | tr -d '\n' | wc -c
}

# Check if the reference genome file exists
if [ -f "$REFERENCE_GENOME_FILE" ]; then
    GENOME_SIZE=$(calculate_genome_size $REFERENCE_GENOME_FILE)
    echo "Calculated genome size: $GENOME_SIZE"

    # Check if the BAM file exists
    if [ -f "$BAM_FILE" ]; then
        BASENAME=$(basename $BAM_FILE .bam)

        macs2 callpeak -t $BAM_FILE -f BAM -g $GENOME_SIZE -n $BASENAME --outdir $OUTPUT_DIR
    else
        echo "BAM file $BAM_FILE not found"
    fi
else
    echo "Reference genome file $REFERENCE_GENOME_FILE not found"
fi
