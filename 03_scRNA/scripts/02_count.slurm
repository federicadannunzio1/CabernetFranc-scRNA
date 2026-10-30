#!/bin/bash
#SBATCH --job-name=cr_count
#SBATCH --output=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/count_%j.out
#SBATCH --error=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/count_%j.err
#SBATCH --partition=dss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --mail-user=federica.dannunzio@uniroma1.it
#SBATCH --mail-type=BEGIN,END,FAIL

PROJ=/lustre/home/dcapozzi/projects/scRNA_protoplasts
export PATH=/lustre/home/dcapozzi/software/cellranger-9.0.1:$PATH

cd $PROJ/03_scRNA/01_alignment

cellranger count \
  --id=VvProt_scRNA \
  --transcriptome=$PROJ/01_reference/cellranger_ref \
  --fastqs=$PROJ/00_raw_data/scRNA \
  --sample=243507_ID4071_1-VvProt-scRNA \
  --localcores=16 \
  --localmem=128 \
  --create-bam=true

echo "count completato"
