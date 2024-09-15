#!/bin/bash
module load samtools
module load deeptools

computeMatrix reference-point --referencePoint center --missingDataAsZero -S m1043.bw m1044.bw m1131_added.bw m1377_1378.bw m1379.bw -R m1043vsm1044_m1043only.bed m1043vsm1044_shared.bed m1043vsm1044_m1044only.bed -a 2000 -b 2000 -o m1043vsm1044_in_m1043_m1044_m1131_added_m1377_m1378_m1379_matrix -bs 20 -p 12
plotHeatmap -m m1043vsm1044_in_m1043_m1044_m1131_added_m1377_m1378_m1379_matrix -o m1043vsm1044_in_m1043_m1044_m1131_added_m1377_m1378_m1379_matrix_heatmap.pdf --boxAroundHeatmaps no --dpi 300 --colorMap Blues --yMin 0 0 0 0 0 --yMax 150 150 100 100 150 --zMin 0 0 0 0 0 --zMax 150 150 100 100 150 --sortRegions descend --sortUsing mean --sortUsingSamples 1 --averageTypeSummaryPlot mean --heatmapHeight 14 --legendLocation upper-left --refPointLabel "center" --xAxisLabel "" --samplesLabel FOXA2_m1043 NKX2-1_m1044 H3k27ac_m1131_added H3k4me1_m1377_1378 H3k4me3_m1379 --regionsLabel "FOXA2 only n=8,716" "shared n=21,817" "NKX2-1 only n=35,181" 
