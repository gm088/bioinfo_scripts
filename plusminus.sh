#!/bin/bash

fname=$1
name=${fname%.bed}
name=${name##*/}

awk '{if($6=="+"){print $0}}' ${fname} > tmp1
awk '{if($6=="-"){print $0}}' ${fname} > tmp2

sort -k1,1 -k2,2n tmp1 > ${name}_plus.bed
sort -k1,1 -k2,2n tmp2 > ${name}_minus.bed

rm tmp1 tmp2
