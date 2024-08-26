#!/bin/bash

module load bedtools

l="../../../loops"


# neither anchors overlap, 
bedtools pairtopair -type neither -a ../w4_fx_nx_linked.bedpe -b $l/w0_mustache_10kb.bedpe > w4_fx_nx_no_w0_loops.bedpe


