#!/bin/bash

if [ -z $3 ]
then
  echo "usage: bedflip.sh fname throwdist offset"
  exit -1
fi

if [ ! -f $1 ]
then
  echo "this file doesn't exist"
  exit -1
fi

echo "okay, make sure you've checked for alternative TSSes"

fname=$1
name=${fname%.bed}
name=${name##*/}
throw=$2
offset=$3

awk -v var="$throw" '{if($6=="+"){print $1"\t"$2-var"\t"$2"\t"$4"_AS""\t"$5"\t""-"}}' ${fname} > int_plus.bed

awk -v var="$throw" '{if($6=="-"){print $1"\t"$3"\t"$3+var"\t"$4"_AS""\t"$5"\t""+"}}' ${fname} > int_minus.bed

cat int_minus.bed int_plus.bed | sort -k1,1 -k2,2n > ${name}_flipped.bed
rm int_minus.bed int_plus.bed

