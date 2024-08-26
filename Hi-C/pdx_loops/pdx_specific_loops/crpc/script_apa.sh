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

h="/projects/b1126/Viriya/hic/lucap/"
dir="plots"
ex="/projects/p20023/Viriya/analysis/foxa2/hi-c/cooltools/output"
if [ ! -d $dir ]; then
	mkdir $dir
fi

for i in CRPC_loop.bedpe
do
j="$(echo "$i" | cut -d'_' -f1,2)"
#j="$(basename "${i}"| sed 's/.[^.]*$//')"
echo ${j}

coolpup.py $h/m1398/LuCaP35CR-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_LuCaP35CR.clpy --features_format bedpe --view $h/cooltools/expected/hg38_arms.bed --expected $h/cooltools/expected/expected_10kb_lucap35CR.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_LuCaP35CR.clpy --output $dir/${j}_LuCaP35CR_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/m1399/LuCaP93-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_LuCaP93.clpy --features_format bedpe --view $h/cooltools/expected/hg38_arms.bed --expected $h/cooltools/expected/expected_10kb_lucap93.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_LuCaP93.clpy --output $dir/${j}_LuCaP93_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/m1400/LuCaP145.2-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_LuCaP145.2.clpy --features_format bedpe --view $h/cooltools/expected/hg38_arms.bed --expected $h/cooltools/expected/expected_10kb_lucap145.2.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_LuCaP145.2.clpy --output $dir/${j}_LuCaP145.2_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/m1401/LuCaP147-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_LuCaP147.clpy --features_format bedpe --view $h/cooltools/expected/hg38_arms.bed --expected $h/cooltools/expected/expected_10kb_lucap147.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_LuCaP147.clpy --output $dir/${j}_LuCaP147_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/m1494/LuCaP70CR-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_LuCaP70CR.clpy --features_format bedpe --view $h/cooltools/expected/hg38_arms.bed --expected $h/cooltools/expected/expected_10kb_lucap70CR.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_LuCaP70CR.clpy --output $dir/${j}_LuCaP70CR_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/m1528/LuCaP77CR-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_LuCaP77CR.clpy --features_format bedpe --view $h/cooltools/expected/hg38_arms.bed --expected $h/cooltools/expected/expected_10kb_lucap77CR.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_LuCaP77CR.clpy --output $dir/${j}_LuCaP77CR_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/m1530/LuCaP145.1-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_LuCaP145.1.clpy --features_format bedpe --view $h/cooltools/expected/hg38_arms.bed --expected $h/cooltools/expected/expected_10kb_lucap145.1.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_LuCaP145.1.clpy --output $dir/${j}_LuCaP145.1_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/m1402/NCI-H660-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_NCI-H660.clpy --features_format bedpe --view $h/cooltools/expected/hg38_arms.bed --expected $h/cooltools/expected/expected_10kb_nci.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_NCI-H660.clpy --output $dir/${j}_NCI-H660_pileup.png --height 3 --vmax 6 --vmin 1

# Timecourse


coolpup.py $h/../lncap/LNCaP-WT-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_W0.clpy --features_format bedpe --view $ex/hg38_arms_lncap.bed --expected $ex/expected_lncap_10k.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_W0.clpy --output $dir/${j}_W0_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/../w2/W2-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_W2.clpy --features_format bedpe --view $ex/hg38_arms.bed --expected $ex/expected_w2_10k.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_W2.clpy --output $dir/${j}_W2_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/../w3/W3-Arima-R1-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_W3.clpy --features_format bedpe --view $ex/hg38_arms.bed --expected $ex/expected_w3_10k.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_W3.clpy --output $dir/${j}_W3_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/../w4/workspace/W4-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_W4.clpy --features_format bedpe --view $ex/hg38_arms.bed --expected $ex/expected_w4_10k.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_W4.clpy --output $dir/${j}_W4_pileup.png --height 3 --vmax 6 --vmin 1

done

export QT_QPA_PLATFORM=offscreen
python stitch_apa.py plots/ CRPC_loop --order LuCaP93 LuCaP145.1 LuCaP145.2 NCI-H660 LuCaP35CR LuCaP70CR LuCaP77CR LuCaP147 W0 W2 W3 W4
