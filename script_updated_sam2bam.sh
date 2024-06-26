#!/bin/bash

#SBATCH --job-name=sam_to_bam     # Job name
#SBATCH --output=sam_to_bam.out   # Standard output and error log
#SBATCH --error=sam_to_bam.err    # Error log
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=4         # Number of CPU cores per task
#SBATCH --mem=4G                  # Total memory for the job
#SBATCH --time=01:00:00           # Time limit hrs:min:sec

# Load the necessary module
module load samtools

# Define the input and output file names
INPUT_SAM="sparus_aurata/SpaAur/SpaAur_ERR11809698_1.fastq.gz.sam"
OUTPUT_BAM="output.bam"
TEMP_SAM="temp.sam"

# Function to filter out problematic lines
filter_sam_file() {
    awk 'BEGIN {OFS="\t"} {
        if ($0 ~ /^@/ || length($10) == length($11)) {
            print $0
        } else {
            print "Filtered out line with mismatched SEQ and QUAL: " NR > "/dev/stderr"
        }
    }' $1 > $2
}

# Validate the input SAM file
samtools quickcheck -v $INPUT_SAM
if [ $? -ne 0 ]; then
  echo "Input SAM file validation failed. Exiting."
  exit 1
fi

# Filter the SAM file
filter_sam_file $INPUT_SAM $TEMP_SAM

# Check if filtering succeeded
if [ $? -ne 0 ]; then
  echo "Error filtering SAM file. Exiting."
  exit 1
fi

# Attempt to convert SAM to BAM
srun samtools view -bS $TEMP_SAM > $OUTPUT_BAM
if [ $? -ne 0 ]; then
  echo "Error converting SAM to BAM. Exiting."
  exit 1
fi

echo "SAM to BAM conversion completed successfully."
