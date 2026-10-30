#!/bin/bash
#SBATCH --job-name=sc_filter
#SBATCH --output=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/filter_%j.out
#SBATCH --error=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/filter_%j.err
#SBATCH --partition=dss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --mail-user=federica.dannunzio@uniroma1.it
#SBATCH --mail-type=BEGIN,END,FAIL

PROJ=/lustre/home/dcapozzi/projects/scRNA_protoplasts

# Usa R dal conda env (NON module load R)
eval "$(conda shell.bash hook)"
conda activate seurat_env

Rscript $PROJ/03_scRNA/scripts/04_filter.R
