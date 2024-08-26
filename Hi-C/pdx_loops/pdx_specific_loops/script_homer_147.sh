#!/bin/bash

module load homer/4.10

sample="147"
res="0.001"

mkdir ${sample}
for i in 93 1451 1452 NCI
do
cp ../${i}v${sample}_${res}/${i}v${sample}.diffloop2.bedpe ${sample}/ 
done


cd ${sample}


merge2Dbed.pl 93v${sample}.diffloop2.bedpe 1451v${sample}.diffloop2.bedpe 1452v${sample}.diffloop2.bedpe NCIv${sample}.diffloop2.bedpe \
-res 10000 -loop -prefix ${sample}_gain > output_${sample}_diffloop2.bedpe 2>&1


