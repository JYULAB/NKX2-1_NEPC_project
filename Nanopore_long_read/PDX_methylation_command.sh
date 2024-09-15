#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00 

module load bowtie2
module load samtools
module load bedtools
module load dimelo

cd /data1/czhao/xiaodong/NCI_H660_DiMeLo/result/pass
dimelo-plot-enrichment-profile -f NCI_H660_DiMeLo_sorted.bam -s NE_CREs AD_CREs -b /data1/czhao/xiaodong/NE_hg38.bed /data1/czhao/xiaodong/AD_hg38.bed -m A -o ./NCI_H660_DiMeLo_NE_AD_Aresult -d 0.05
dimelo-plot-enrichment-profile -f NCI_H660_DiMeLo_sorted.bam -s NE_CREs AD_CREs -b /data1/czhao/xiaodong/NE_hg38.bed /data1/czhao/xiaodong/AD_hg38.bed -m CG -o ./NCI_H660_DiMeLo_NE_AD_CGresult -d 0.05
