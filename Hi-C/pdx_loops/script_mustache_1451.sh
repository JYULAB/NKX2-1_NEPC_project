#!/bin/bash
#SBATCH -A p20023
#SBATCH -p normal 
#SBATCH -n 1
#SBATCH --mem=40G
#SBATCH -t 48:00:00
#SBATCH --mail-user=viriyakeo2023@u.northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=z_mustache_93
#SBATCH -o "%x.o%j"

source ~/.bashrc
conda activate mustache
pt2="0.001"
samples=("35CR" "147" "70CR" "77CR")
num=("m1398" "m1401" "m1494" "m1528")
main="1451"

for ((index=0;index<4;index++)) 
do
j="${main}v${samples[$index]}"
mkdir ${j}_${pt2}

echo "index ${index}"
echo "j ${j}"
echo "folder ${j}_${pt2}"
echo "CRPC ${num[$index]}/LuCaP${samples[$index]}-Arima-allReps-filtered.mcool"

/home/vkp2256/mustache/mustache/diff_mustache.py -f1 ../../m1530/LuCaP145.1-Arima-allReps-filtered.mcool \
-f2 ../../${num[$index]}/LuCaP${samples[$index]}-Arima-allReps-filtered.mcool -pt 0.05 -pt2 $pt2 \
-o ${j}_${pt2}/${j} \
-r 10000 -st 0.8



for i in ${j}_${pt2}/${j}.diffloop1 ${j}_${pt2}/${j}.diffloop2 ${j}_${pt2}/${j}.loop1 ${j}_${pt2}/${j}.loop2
do
tail -n+2 ${i} > ${j}_${pt2}/${i}.bedpe
done

wc -l ${j}_${pt2}/*bedpe > ${j}_${pt2}/mustache_counts.txt

done
conda deactivate
