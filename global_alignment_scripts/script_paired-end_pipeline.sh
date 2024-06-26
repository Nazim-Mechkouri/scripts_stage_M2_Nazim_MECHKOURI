#!/bin/bash
#SBATCH --job-name=atacseq_pipeline
#SBATCH --output=atacseq_pipeline_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# Load modules
module load trimmomatic
module load bowtie2
module load samtools

# Input parameters
READ1=$1
READ2=$2
REFERENCE_GENOME=$3
ADAPTERS=$4
OUTPUT_DIR=$5
OUTPUT_PREFIX=$6

########################################################################################
# Step 0 : pre-processing

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Define the Bowtie 2 index prefix
INDEX_PREFIX=${REFERENCE_GENOME%.*}_index

# Check that the Bowtie 2 index exists and we are having 6 index files.
#if [ ! -f ${INDEX_PREFIX}.1.bt2 ] || [ ! -f ${INDEX_PREFIX}.2.bt2 ] || \
#   [ ! -f ${INDEX_PREFIX}.3.bt2 ] || [ ! -f ${INDEX_PREFIX}.4.bt2 ] || \
#   [ ! -f ${INDEX_PREFIX}.rev.1.bt2 ] || [ ! -f ${INDEX_PREFIX}.rev.2.bt2 ]; then
#  echo "Error: Bowtie 2 index not found. Please build the index first."
#  exit 1
#fi


# Extract the base name of the file for output naming
BASE_NAME=$(basename $READ1 _R1.fastq)
FILE_OUTPUT_PREFIX=$OUTPUT_DIR/${OUTPUT_PREFIX}_${BASE_NAME}


#######################################################################################

# Step 1: Trimming via Trimmomatics (mostly used in biblio)
echo "Starting trimming for $READ1 and $READ2..."
trimmomatic PE -threads 8 \
    $READ1 $READ2 \
    ${FILE_OUTPUT_PREFIX}_trimmed_R1_paired.fastq ${FILE_OUTPUT_PREFIX}_trimmed_R1_unpaired.fastq \
    ${FILE_OUTPUT_PREFIX}_trimmed_R2_paired.fastq ${FILE_OUTPUT_PREFIX}_trimmed_R2_unpaired.fastq \
    ILLUMINACLIP:$ADAPTERS:2:30:10 SLIDINGWINDOW:4:20 MINLEN:25
echo "Trimming completed for $READ1 and $READ2."


#######################################################################################
# Step 2: Alignment with Bowtie2 : pipeline ofr paired-end data
echo "Starting alignment for ${FILE_OUTPUT_PREFIX}_trimmed_R1_paired.fastq and ${FILE_OUTPUT_PREFIX}_trimmed_R2_paired.fastq..."
bowtie2 -p 8 -x $INDEX_PREFIX \
    -1 ${FILE_OUTPUT_PREFIX}_trimmed_R1_paired.fastq \
    -2 ${FILE_OUTPUT_PREFIX}_trimmed_R2_paired.fastq \
    -S ${FILE_OUTPUT_PREFIX}.sam
echo "Alignment completed for ${FILE_OUTPUT_PREFIX}_trimmed_R1_paired.fastq and ${FILE_OUTPUT_PREFIX}_trimmed_R2_paired.fastq."


######################### SAMTools and conversions ####################################
# Step 3: Convert SAM to BAM
echo "Converting SAM to BAM for $FILE_OUTPUT_PREFIX..."
samtools view -bS ${FILE_OUTPUT_PREFIX}.sam > ${FILE_OUTPUT_PREFIX}.bam
echo "Conversion to BAM completed for $FILE_OUTPUT_PREFIX."

# Step 4: Sort BAM
echo "Sorting BAM for $FILE_OUTPUT_PREFIX..."
samtools sort ${FILE_OUTPUT_PREFIX}.bam -o ${FILE_OUTPUT_PREFIX}_sorted.bam
echo "Sorting completed for $FILE_OUTPUT_PREFIX."

# Step 5: Index BAM
echo "Indexing BAM for $FILE_OUTPUT_PREFIX..."
samtools index ${FILE_OUTPUT_PREFIX}_sorted.bam
echo "Indexing completed for $FILE_OUTPUT_PREFIX."

# Cleanup intermediate files
rm ${FILE_OUTPUT_PREFIX}.sam ${FILE_OUTPUT_PREFIX}.bam

echo "Pipeline completed successfully for $READ1 and $READ2."
