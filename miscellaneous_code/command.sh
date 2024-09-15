PATHTO=/data1/czhao/xiaodong/FOXA2/ROSE
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin
ROSE_main.py -g HG19 -i m1407_m1390asCtrl_peaks.narrowPeak.gff -r m1407_sorted.bam -o m1407_m1390asCtrl_SE
