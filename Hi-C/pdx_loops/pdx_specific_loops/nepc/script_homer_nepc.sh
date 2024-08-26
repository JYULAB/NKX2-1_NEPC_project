#!/bin/bash

module load homer/4.10
mkdir results

for i in *bedpe
do

j=${i%%_*}

cp ${i} results/${j}.bedpe
done

cd results

merge2Dbed.pl 93.bedpe 1451.bedpe 1452.bedpe nci.bedpe \
-res 10000 -loop -prefix nepc_gain > output_nepc.bedpe 2>&1





