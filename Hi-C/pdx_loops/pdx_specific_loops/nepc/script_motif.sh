#!/bin/bash

#SBATCH -A b1042 
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mem=30G
#SBATCH --mail-user=viriyakeo2023@u.northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=z_motif_atac_anchor
#SBATCH -o "%x.o%j" 

module purge

# load modules
module load python/anaconda3.6
module load homer


for i in *bed
do

echo $i
findMotifsGenome.pl ${i} /projects/p20023/Viriya/software/Homer4.10/data/genomes/hg38/ ${i%.*} -size 200

done

