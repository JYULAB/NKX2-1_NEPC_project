#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00 

module load bowtie2
module load samtools
module load bedtools
module load dimelo


awk -F"\t" '$11!="nan"&& $11>0' RRMS_FOXA2_2day_IgG_5mc.bed > RRMS_FOXA2_2day_IgG_5mc_remove_nan_zero.bed
awk -v OFS='\t' '{print $1,$2,$3,$11}' RRMS_FOXA2_2day_IgG_5mc_remove_nan_zero.bed | sort -k1,1 -k2,2n > RRMS_FOXA2_2day_IgG_5mc_remove_nan_zero.bedgraph
/projects/b1042/YuLab/./bedGraphToBigWig RRMS_FOXA2_2day_IgG_5mc_remove_nan_zero.bedgraph /projects/b1042/YuLab/nanopore/reference_genome/hg38/hg38.chrom.sizes RRMS_FOXA2_2day_IgG_5mc_remove_nan_zero.bigwig


dm.plot_browser(bam, sampleName, "chr11:2086423-2091187", "A+CG", outDir, threshA=153, threshC=153, static=True, smooth=100, min_periods=10)
