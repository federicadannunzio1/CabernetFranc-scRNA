#!/bin/bash
#SBATCH --job-name=cr_mkref
#SBATCH --output=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/mkref_%j.out
#SBATCH --error=/lustre/home/dcapozzi/projects/scRNA_protoplasts/logs/mkref_%j.err
#SBATCH --partition=dss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --mail-user=federica.dannunzio@uniroma1.it
#SBATCH --mail-type=BEGIN,END,FAIL

PROJ=/lustre/home/dcapozzi/projects/scRNA_protoplasts
export PATH=/lustre/home/dcapozzi/software/cellranger-9.0.1:$PATH

cd $PROJ/01_reference

# Filtra il GTF per Cell Ranger
cellranger mkgtf \
  $PROJ/01_reference/PN40024_5.1_on_T2T_ref_with_names_no_chr00_fixed.gtf \
  $PROJ/01_reference/PN40024_5.1_on_T2T_ref_cellranger.gtf

# Genera il reference
cellranger mkref \
  --genome=vitis_T2T_no_chr00 \
  --fasta=$PROJ/01_reference/T2T_ref_no_chr00.fasta \
  --genes=$PROJ/01_reference/PN40024_5.1_on_T2T_ref_cellranger.gtf \
  --nthreads=16 \
  --memgb=64

# Sposta nella posizione finale
mv $PROJ/01_reference/mkref_vitis_T2T_no_chr00/outs/vitis_T2T_no_chr00 $PROJ/01_reference/cellranger_ref

echo "mkref completato"
