#!/bin/bash
#SBATCH --job-name=sc_markers
#SBATCH --output=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/markers_%j.out
#SBATCH --error=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/markers_%j.err
#SBATCH --partition=dss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=03:00:00
#SBATCH --mail-user=federica.dannunzio@uniroma1.it
#SBATCH --mail-type=BEGIN,END,FAIL

PROJ=/lustre/home/dcapozzi/projects/scRNA_protoplasts

eval "$(conda shell.bash hook)"
conda activate seurat_env

Rscript $PROJ/repo/03_scRNA/scripts/06_markers.R
