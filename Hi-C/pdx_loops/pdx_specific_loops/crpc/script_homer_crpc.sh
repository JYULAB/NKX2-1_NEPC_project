#!/bin/bash

module load homer/4.10
mkdir results

for i in *bedpe
do

j=${i%%_*}

cp ${i} results/${j}.bedpe
done

cd results

merge2Dbed.pl 35CR.bedpe 70CR.bedpe 77CR.bedpe 147.bedpe \
-res 10000 -loop -prefix crpc_gain > output_crpc.bedpe 2>&1





