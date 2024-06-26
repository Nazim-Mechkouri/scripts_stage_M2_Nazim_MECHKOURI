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
BAM_PARENT_DIR=$1
REFERENCE_GENOME_PARENT_DIR=$2
OUTPUT_DIR=$3

# Create my directory if not existing already
mkdir -p $OUTPUT_DIR

# Get genome size from reference_genome, using only wc -c and not counting the header. The genome size is neededfor MACS2
calculate_genome_size() {
    local genome_file=$1
    grep -v "^>" $genome_file | tr -d '\n' | wc -c
}

# Peak calling for each BAM file in each subdirectory
for SPECIES_DIR in $BAM_PARENT_DIR/*/; do
    SPECIES_NAME=$(basename $SPECIES_DIR)
    REFERENCE_GENOME_FILE=$(find $REFERENCE_GENOME_PARENT_DIR/$SPECIES_NAME -name "*.fna" | head -n 1)

    if [ -f "$REFERENCE_GENOME_FILE" ]; then
        GENOME_SIZE=$(calculate_genome_size $REFERENCE_GENOME_FILE)
        echo "Calculated genome size for $SPECIES_NAME: $GENOME_SIZE"

        for BAM_FILE in $SPECIES_DIR/*.bam; do
            BASENAME=$(basename $BAM_FILE .bam)
            SPECIES_OUTPUT_DIR=$OUTPUT_DIR/$SPECIES_NAME

            mkdir -p $SPECIES_OUTPUT_DIR

            macs2 callpeak -t $BAM_FILE -f BAM -g $GENOME_SIZE -n $BASENAME --outdir $SPECIES_OUTPUT_DIR
        done
    else
        echo "Reference genome for $SPECIES_NAME not found in $REFERENCE_GENOME_PARENT_DIR/$SPECIES_NAME"
    fi
done

echo "Peak calling completed !"
