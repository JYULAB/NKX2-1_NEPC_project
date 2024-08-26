#!/bin/bash

#SBATCH -A p20023
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH --mem=20G
#SBATCH --mail-user=viriyakeo2023@u.northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=z_apa
#SBATCH -o "%x.o%j"

source ~/.bashrc
conda activate peakachu

h="/projects/p20023/Viriya/analysis/foxa2/hi-c"
dir="flanks"
if [ ! -d $dir ]; then
	mkdir $dir
fi

for i in *loops.bedpe
do
echo ${i}

j="$(basename "${i}"| sed 's/.[^.]*$//')"

#coolpup.py $h/mcool/W4-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_w4.clpy --features_format bedpe --view $h/cooltools/output/hg38_arms.bed --expected $h/cooltools/output/expected_w4_10k.tsv --flank 80000

#plotpup.py --input_pups $dir/${j}_w4.clpy --output $dir/${j}_w4_pileup.png --height 3 --vmax 6
#plotpup.py --input_pups $dir/${j}_w4.clpy --output $dir/${j}_w4_pileup_1.png --height 3 --vmax 6 --vmin 1
plotpup.py --input_pups $dir/${j}_w4.clpy --output $dir/${j}_w4_pileup.pdf --height 3 --vmax 6 --vmin 1
#plotpup.py --input_pups $dir/${j}_w4.clpy --output $dir/${j}_w4_pileup_2.png --height 3 --vmax 5 --vmin 0.5

#coolpup.py $h/mcool/W3-Arima-R1-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_w3.clpy --features_format bedpe --view $h/cooltools/output/hg38_arms.bed --expected $h/cooltools/output/expected_w3_10k.tsv --flank 80000
#plotpup.py --input_pups $dir/${j}_w3.clpy --output $dir/${j}_w3_pileup.png --height 3 --vmax 6
#plotpup.py --input_pups $dir/${j}_w3.clpy --output $dir/${j}_w3_pileup_1.png --height 3 --vmax 6 --vmin 1
plotpup.py --input_pups $dir/${j}_w3.clpy --output $dir/${j}_w3_pileup.pdf --height 3 --vmax 6 --vmin 1
#plotpup.py --input_pups $dir/${j}_w3.clpy --output $dir/${j}_w3_pileup_2.png --height 3 --vmax 5 --vmin 0.5

#coolpup.py $h/mcool/W2-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_w2.clpy --features_format bedpe --view $h/cooltools/output/hg38_arms.bed --expected $h/cooltools/output/expected_w2_10k.tsv --flank 80000
#plotpup.py --input_pups $dir/${j}_w2.clpy --output $dir/${j}_w2_pileup.png --height 3 --vmax 6
#plotpup.py --input_pups $dir/${j}_w2.clpy --output $dir/${j}_w2_pileup_1.png --height 3 --vmax 6 --vmin 1
plotpup.py --input_pups $dir/${j}_w2.clpy --output $dir/${j}_w2_pileup.pdf --height 3 --vmax 6 --vmin 1
#plotpup.py --input_pups $dir/${j}_w2.clpy --output $dir/${j}_w2_pileup_2.png --height 3 --vmax 5 --vmin 0.5

#coolpup.py $h/mcool/LNCaP-WT-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}.clpy --features_format bedpe --view $h/cooltools/output/hg38_arms_lncap.bed --expected $h/cooltools/output/expected_lncap_10k.tsv --flank 80000
#plotpup.py --input_pups $dir/${j}.clpy --output $dir/${j}_lncap_pileup.png --height 3 --vmax 6
#plotpup.py --input_pups $dir/${j}.clpy --output $dir/${j}_lncap_pileup_1.png --height 3 --vmax 6 --vmin 1
plotpup.py --input_pups $dir/${j}.clpy --output $dir/${j}_lncap_pileup.pdf --height 3 --vmax 6 --vmin 1
#plotpup.py --input_pups $dir/${j}.clpy --output $dir/${j}_lncap_pileup_2.png --height 3 --vmax 5 --vmin 0.5



done

