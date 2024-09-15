#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00 

module load bowtie2
module load samtools
module load bedtools
module load dimelo

cd /data1/czhao/xiaodong/FOXA2_D2_H3k4me1_dual/result/pass
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q4 q3 q2 q1 -b FOXA2_ChIP_peak_4w_only_peak_q4.bed FOXA2_ChIP_peak_4w_only_peak_q3.bed FOXA2_ChIP_peak_4w_only_peak_q2.bed FOXA2_ChIP_peak_4w_only_peak_q1.bed -m A -o ./FOXA2_D2_H3k4me1_dual_Aresult -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q4 q3 q2 q1 -b FOXA2_ChIP_peak_4w_only_peak_q4.bed FOXA2_ChIP_peak_4w_only_peak_q3.bed FOXA2_ChIP_peak_4w_only_peak_q2.bed FOXA2_ChIP_peak_4w_only_peak_q1.bed -m CG -o ./FOXA2_D2_H3k4me1_dual_CGresult -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q1 -b FOXA2_ChIP_peak_4w_only_peak_q1.bed -m A -o ./FOXA2_D2_H3k4me1_dual_Aresult_q1 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q2 -b FOXA2_ChIP_peak_4w_only_peak_q2.bed -m A -o ./FOXA2_D2_H3k4me1_dual_Aresult_q2 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q3 -b FOXA2_ChIP_peak_4w_only_peak_q3.bed -m A -o ./FOXA2_D2_H3k4me1_dual_Aresult_q3 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q4 -b FOXA2_ChIP_peak_4w_only_peak_q4.bed -m A -o ./FOXA2_D2_H3k4me1_dual_Aresult_q4 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q1 -b FOXA2_ChIP_peak_4w_only_peak_q1.bed -m CG -o ./FOXA2_D2_H3k4me1_dual_CGresult_q1 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q2 -b FOXA2_ChIP_peak_4w_only_peak_q2.bed -m CG -o ./FOXA2_D2_H3k4me1_dual_CGresult_q2 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q3 -b FOXA2_ChIP_peak_4w_only_peak_q3.bed -m CG -o ./FOXA2_D2_H3k4me1_dual_CGresult_q3 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D2_H3k4me1_dual_sorted.bam -s q4 -b FOXA2_ChIP_peak_4w_only_peak_q4.bed -m CG -o ./FOXA2_D2_H3k4me1_dual_CGresult_q4 -d 0.05

cd /data1/czhao/xiaodong/FOXA2_D28_H3k4me1_dual/result/pass
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q4 q3 q2 q1 -b FOXA2_ChIP_peak_4w_only_peak_q4.bed FOXA2_ChIP_peak_4w_only_peak_q3.bed FOXA2_ChIP_peak_4w_only_peak_q2.bed FOXA2_ChIP_peak_4w_only_peak_q1.bed -m A -o ./FOXA2_D28_H3k4me1_dual_Aresult -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q4 q3 q2 q1 -b FOXA2_ChIP_peak_4w_only_peak_q4.bed FOXA2_ChIP_peak_4w_only_peak_q3.bed FOXA2_ChIP_peak_4w_only_peak_q2.bed FOXA2_ChIP_peak_4w_only_peak_q1.bed -m CG -o ./FOXA2_D28_H3k4me1_dual_CGresult -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q1 -b FOXA2_ChIP_peak_4w_only_peak_q1.bed -m A -o ./FOXA2_D28_H3k4me1_dual_Aresult_q1 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q2 -b FOXA2_ChIP_peak_4w_only_peak_q2.bed -m A -o ./FOXA2_D28_H3k4me1_dual_Aresult_q2 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q3 -b FOXA2_ChIP_peak_4w_only_peak_q3.bed -m A -o ./FOXA2_D28_H3k4me1_dual_Aresult_q3 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q4 -b FOXA2_ChIP_peak_4w_only_peak_q4.bed -m A -o ./FOXA2_D28_H3k4me1_dual_Aresult_q4 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q1 -b FOXA2_ChIP_peak_4w_only_peak_q1.bed -m CG -o ./FOXA2_D28_H3k4me1_dual_CGresult_q1 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q2 -b FOXA2_ChIP_peak_4w_only_peak_q2.bed -m CG -o ./FOXA2_D28_H3k4me1_dual_CGresult_q2 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q3 -b FOXA2_ChIP_peak_4w_only_peak_q3.bed -m CG -o ./FOXA2_D28_H3k4me1_dual_CGresult_q3 -d 0.05
dimelo-plot-enrichment-profile -f FOXA2_D28_H3k4me1_dual_sorted.bam -s q4 -b FOXA2_ChIP_peak_4w_only_peak_q4.bed -m CG -o ./FOXA2_D28_H3k4me1_dual_CGresult_q4 -d 0.05
