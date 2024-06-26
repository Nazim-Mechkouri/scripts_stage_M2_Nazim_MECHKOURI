#!/bin/bash

#SBATCH --job-name=MafHomoCanis
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem 40GB                     # mémoire vive pour l'ensemble des cœurs
#SBATCH -o /shared/projects/pten_enh_phylo/project_pten_nazim/synteny/slurm_output/slurm.%N.%j.out           # STDOUT
#SBATCH -e /shared/projects/pten_enh_phylo/project_pten_nazim/synteny/slurm_output/slurm.%N.%j.err           # STDERR



################################################################################
# This scripts filters mulputiple alignment files .maf                         #
#                                                                              #
# If problems encoutered concerning the input or usage of the script           #
# please refer to the help section below or use : ./mafReader -h               #
#                                                                              #
#  version updates and perspectives :-                                         #
#                                    -                                         #
#                                    -                                         #
#                                    -                                         #
#                                    -                                         #
################################################################################

################################################################################
# Help function                                                                #
################################################################################
Help()
{
   # Display Help
   RED='\033[1;33m'
   NC='\033[0m'

   echo
   echo "####	This script takes 6 arguments as follows	####"
   echo
   echo
   echo -e "${RED}Usage : ${NC}./mafReader.sh MAFfile SpecieGenome ChrNumber Start-Coordinate End-Coordinate OutputFile"
   echo
   echo -e "${RED}-mafFile :${NC} obtained from UCSC genome browser data"
   echo -e "${RED}-SpecieGenome & ChrNumber :${NC} name of the genome.assembly as it is referenced in the maf file (ex : hg38,mmus19,CanFam5,BosTau9 ...etc). ChrNumber only ("10" and not "chr10")"
   echo -e "${RED}-Start-Coordinate :${NC} Start of the targeted sequence to be looked for"
   echo -e "${RED}-End-Coordinate :${NC} End of the targeted sequence to be looked for"
}


################################################################################
# Process the input options. Add options as needed.                            #
################################################################################

# Get the options
while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
     \?) # incorrect option
         echo "Error: Invalid option"
         exit;;
   esac
done


################################################################################
#                              Main program                                    #
################################################################################

# Define output file
output_file=$6

# Initialize the output file
> "$output_file"


# my inputs #
maf_file=$1
input_specie=$2
input_chrom=$3
sequence_start_position=$4
sequence_end_position=$5

# my counts #
found_specie=0
found_chrom=0
found_candidate_sequence=0
cpt=0
temporary_count=1001

# split the filename from the Path #
p=$maf_file
file=$(basename "$p")

# Initialization and quality of life
echo "The current file is: $file"
echo "Species of interest: $input_specie"
echo "Chromosome of interest: chr$input_chrom"
echo "Looking near position: $sequence_start_position and/or $sequence_end_position"
echo
echo "Searching for a match between the current position and the given position in the file..."

###################################################################################################
# Reading the file
###################################################################################################

# Flag variable to signal when required data has been found
found_required_data=0  

while IFS= read -r line || [[ -n $line ]]; do
    # Check if the line starts with 'a' indicating a NEW alignment block #
    if [[ $line == "a"* ]]; then
        found_specie=0
        found_chrom=0
        found_candidate_sequence=0
        block_info="$line\n"  # Initialize block_info variable with the line
    fi

    # Check if the line starts with 's' indicating a sequence line
    if [[ $line == "s"* ]]; then
        # Extract the sequence from the line using an array list #
        line_parts=($line)
        current_specie=$(echo ${line_parts[1]} | cut -d "." -f 1)  # specie of the current line (for pairwise human/dog for example its either hg38 or CanFam)
        current_chrom=$(echo ${line_parts[1]} | cut -d "." -f 2)  # chromosome where the current sequence is located
        current_start_position=$(echo ${line_parts[2]} | cut -f 2) 

        cpt=$((cpt+1))
        if [[ $cpt == $temporary_count ]]; then 
            echo $cpt $current_chrom $current_specie $current_start_position
            temporary_count=$((temporary_count*2-1))
        fi    
        
        if [[ $current_start_position -ge $sequence_start_position && $current_start_position -le $sequence_end_position ]]; then
            found_candidate_sequence=1

            if [[ $current_specie == $input_specie && $current_chrom == "chr$input_chrom" ]]; then
                found_specie=1
                found_chrom=1
                echo "Found a match for the given sequence coordinates:"
                echo "Current position: $current_start_position matches with the given input position: $sequence_start_position"
            fi
        fi    
    fi

    # Check if all conditions are met, then append block_info to the output file
    if [[ $found_specie -eq 1 && $found_chrom -eq 1 && $found_candidate_sequence -eq 1 ]]; then
        echo -e "$block_info\n$line" >> "$output_file"
        found_required_data=1  # Set flag to 1 once required data is found
    fi

    # Append block info to block_info variable
    block_info="$block_info\n$line"

done < "$maf_file"  

# Check if required data was found
if [[ $found_required_data -eq 0 ]]; then
    echo "Required data not found in the input file."
    exit 1  # Or any other appropriate action
fi

echo "Script execution completed."
