#!/bin/bash

#SBATCH --job-name=HomoMus_align
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem 40GB                     # mémoire vive pour l'ensemble des cœurs
#SBATCH -o /shared/projects/pten_enh_phylo/project_pten_nazim/synteny/slurm_output/slurm.%N.%j.out           # STDOUT
#SBATCH -e /shared/projects/pten_enh_phylo/project_pten_nazim/synteny/slurm_output/slurm.%N.%j.err           # STDERR

module load minimap2/2.24
perl ../NGenomeSyn/bin/NGenomeSyn  -InConf ../slurm_output/final_outputs/run_homo_mus_1/homo_gal.conf -OutPut red_orange
