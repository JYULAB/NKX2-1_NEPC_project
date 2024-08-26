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



h="/projects/b1126/Viriya/hic/"
dir="plots"
ex="/projects/p20023/Viriya/analysis/foxa2/hi-c/cooltools/output"
ex_nci="/projects/b1126/Viriya/hic/lucap/cooltools/expected"
if [ ! -d $dir ]; then
	mkdir $dir
fi
sample="nci"
for i in w0_vs_w4.diffloop1.bedpe w0_vs_w4.diffloop2.bedpe
do
echo ${i}
j="$(basename "${i}" | sed 's/\.[^.]*$//')"
echo ${j}

coolpup.py $h/lncap/LNCaP-WT-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_W0.clpy --features_format bedpe --view $ex/hg38_arms_lncap.bed --expected $ex/expected_lncap_10k.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_W0.clpy --output $dir/${j}_W0_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/w2/W2-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_W2.clpy --features_format bedpe --view $ex/hg38_arms.bed --expected $ex/expected_w2_10k.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_W2.clpy --output $dir/${j}_W2_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/w3/W3-Arima-R1-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_W3.clpy --features_format bedpe --view $ex/hg38_arms.bed --expected $ex/expected_w3_10k.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_W3.clpy --output $dir/${j}_W3_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/w4/workspace/W4-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_W4.clpy --features_format bedpe --view $ex/hg38_arms.bed --expected $ex/expected_w4_10k.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_W4.clpy --output $dir/${j}_W4_pileup.png --height 3 --vmax 6 --vmin 1

coolpup.py $h/lucap/m1402/NCI-H660-Arima-allReps-filtered.mcool::resolutions/10000 ${i} --nproc 2 -o $dir/${j}_nci.clpy --features_format bedpe --view $h/lucap/cooltools/expected/hg38_arms.bed --expected $h/lucap/cooltools/expected/expected_10kb_nci.tsv --flank 50000

plotpup.py --input_pups $dir/${j}_nci.clpy --output $dir/${j}_nci_pileup.png --height 3 --vmax 6 --vmin 1
done


python ../../stitch_apa.py plots/ w0_vs_w4.diffloop1 --order W0 W2 W3 W4 nci
python ../../stitch_apa.py plots/ w0_vs_w4.diffloop2 --order W0 W2 W3 W4 nci

conda deactivate
