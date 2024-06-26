#!/bin/bash

#SBATCH --job-name=HomoMus_align
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem 50GB                     # mémoire vive pour l'ensemble des cœurs
#SBATCH -o /shared/projects/pten_enh_phylo/project_pten_nazim/synteny/slurm_output/slurm.%N.%j.out           # STDOUT
#SBATCH -e /shared/projects/pten_enh_phylo/project_pten_nazim/synteny/slurm_output/slurm.%N.%j.err           # STDERR

infile1=$1
infile2=$2
outfile=$3

module load mummer4/4.0.0rc1
module load perl/5.26.2

perl ../../NGenomeSyn/bin/GetTwoGenomeSyn.pl  -InGenomeA $infile1 -InGenomeB $infile2 -BinDir . -MinLenA 1000000 -MinlenB 1000000 -MinAlnLen 2000 -OutPrefix $outfile
