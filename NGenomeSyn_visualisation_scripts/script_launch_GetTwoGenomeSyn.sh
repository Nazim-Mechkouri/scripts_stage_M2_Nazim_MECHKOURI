#!/bin/bash

#SBATCH --job-name=HomoMus_align
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem 40GB                     # mémoire vive pour l'ensemble des cœurs
#SBATCH -o /shared/projects/pten_enh_phylo/project_pten_nazim/synteny/slurm_output/slurm.%N.%j.out           # STDOUT
#SBATCH -e /shared/projects/pten_enh_phylo/project_pten_nazim/synteny/slurm_output/slurm.%N.%j.err           # STDERR

module load minimap2/2.24
perl ../NGenomeSyn/bin/GetTwoGenomeSyn.pl  -InGenomeA ../reference_genomes/homo_sapiens/homo_sapiens_GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna -InGenomeB ../reference_genomes/mus_musculus/mus_musculus_GCA_000001635.9_GRCm39_genomic.fna -MappingBin minimap2 -MappingDir . -OutPrefix homo_gal -MappingPara "-x asm20 -k17"

