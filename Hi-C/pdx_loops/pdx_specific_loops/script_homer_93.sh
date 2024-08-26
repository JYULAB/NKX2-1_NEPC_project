#!/bin/bash

module load homer/4.10

sample="93"
res="0.001"

mkdir ${sample}

cp ../${sample}v35CR_${res}/${sample}v35CR.diffloop1.bedpe ${sample}/${sample}v35CR.diffloop1.bedpe
cp ../${sample}v70CR_${res}/${sample}v70CR.diffloop1.bedpe ${sample}/${sample}v70CR.diffloop1.bedpe
cp ../${sample}v77CR_${res}/${sample}v77CR.diffloop1.bedpe ${sample}/${sample}v77CR.diffloop1.bedpe
cp ../${sample}v147_${res}/${sample}v147.diffloop1.bedpe ${sample}/${sample}v147.diffloop1.bedpe


cp ../${sample}v35CR_${res}/${sample}v35CR.diffloop2.bedpe ${sample}/${sample}v35CR.diffloop2.bedpe 
cp ../${sample}v70CR_${res}/${sample}v70CR.diffloop2.bedpe ${sample}/${sample}v70CR.diffloop2.bedpe
cp ../${sample}v77CR_${res}/${sample}v77CR.diffloop2.bedpe ${sample}/${sample}v77CR.diffloop2.bedpe 
cp ../${sample}v147_${res}/${sample}v147.diffloop2.bedpe ${sample}/${sample}v147.diffloop2.bedpe 


cd ${sample}

merge2Dbed.pl ${sample}v35CR.diffloop1.bedpe ${sample}v70CR.diffloop1.bedpe ${sample}v77CR.diffloop1.bedpe ${sample}v147.diffloop1.bedpe -res 10000 -loop -prefix 93_gain > output_93_diffloop1.bedpe 2>&1

merge2Dbed.pl ${sample}v35CR.diffloop2.bedpe ${sample}v70CR.diffloop2.bedpe ${sample}v77CR.diffloop2.bedpe ${sample}v147.diffloop2.bedpe -res 10000 -loop -prefix 93_loss > output_93_diffloop2.bedpe 2>&1


