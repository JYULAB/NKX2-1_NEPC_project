#!/bin/bash
#SBATCH -A p20023
#SBATCH -p normal  
#SBATCH -n 1
#SBATCH --mem=40G
#SBATCH -t 48:00:00
#SBATCH --mail-user=viriyakeo2023@u.northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=z_mustache
#SBATCH -o "%x.o%j"

source ~/.bashrc
conda activate mustache


/home/vkp2256/mustache/mustache/diff_mustache.py -f1 ../../../../lncap/LNCaP-WT-Arima-allReps-filtered.mcool \
-f2 ../../../../w4/workspace/W4-Arima-allReps-filtered.mcool -pt 0.05 -pt2 0.001 -o w0_vs_w4 \
-r 10000 -st 0.8


conda deactivate
for i in *.diffloop1 *.diffloop2 *.loop1 *.loop2

do
tail -n+2 ${i} > ${i}.bedpe

done

wc -l *bedpe > num_loops.txt


