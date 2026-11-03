#!/bin/bash
#SBATCH --job-name=sc_qc
#SBATCH --output=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/qc_%j.out
#SBATCH --error=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/qc_%j.err
#SBATCH --partition=dss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --mail-user=federica.dannunzio@uniroma1.it
#SBATCH --mail-type=BEGIN,END,FAIL

PROJ=/lustre/home/dcapozzi/projects/scRNA_protoplasts

# Usa R dal conda env (NON module load R — il modulo non vede i pacchetti conda)
eval "$(conda shell.bash hook)"
conda activate seurat_env

Rscript $PROJ/repo/03_scRNA/scripts/03_qc.R
