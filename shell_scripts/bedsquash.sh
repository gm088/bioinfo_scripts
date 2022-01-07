#!/bin/bash

#squash to TSS or TES
if [ -z $2 ]
then
  echo "usage: bedsquash.sh fname TSS/TES/MID"
  exit -1
fi

if [ ! -f $1 ]
then
  echo "this file doesn't exist"
  exit -1
fi

fname=$1
name=${fname%.bed}
name=${name##*/}

if [[ $2 == "TES" ]]
then
	##squash to 3'
	
	awk -F'\t' -vOFS='\t' '{if($6=="+"){$2=$3-1}{print $0}}' $1 \
	| awk -F'\t' -vOFS='\t' '{if($6=="-"){$3=$2+1}{print $0}}' \
	| sort -k1,1 -k2,2n > ${name}_3p_squash.bed
	
fi

if [[ $2 == "TSS" ]]
then
	awk -F'\t' -vOFS='\t' '{if($6=="+"){$3=$2+1}{print $0}}' $1 \
	| awk -F'\t' -vOFS='\t' '{if($6=="-"){$2=$3-1}{print $0}}' \
	| sort -k1,1 -k2,2n > ${name}_5p_squash.bed
	
fi

if [[ $2 == "MID" ]]
then
	##midpoint
	
	awk '{printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1, ($2+$3)/2, (($2+$3)/2)+1, $4, $5, $6}' > ${name}_midpoint.bed
	
fi

