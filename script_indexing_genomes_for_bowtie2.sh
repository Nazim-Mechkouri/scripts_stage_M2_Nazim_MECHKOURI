#!/bin/bash
#SBATCH --job-name=build_bowtie2_index
#SBATCH --output=build_bowtie2_index_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# Load Bowtie2 module
module load bowtie2

# Input parameters
REFERENCE_GENOME=$1
INDEX_PREFIX=${REFERENCE_GENOME%.*}_index

# Build Bowtie2 index. Index must be 6 (the other script checks if this information is true), and are needed to align via Bowtie2
echo "Building Bowtie 2 index for $REFERENCE_GENOME..."
bowtie2-build $REFERENCE_GENOME $INDEX_PREFIX
echo "Indexing completed."
