#!/bin/bash

fname=$1
name=${fname%.bed}
name=${name##*/}


if [ -z $3 ]
then
	echo "usage: bedfile binsize strand"
	exit 1
fi

bed=$1

if [[ $3 == "minus" ]]
then
#  bed="${name}_minus.bed"
  bigwigs=`ls ../../bigwigs_norm/*_rev.bw`
  strand="minus"
elif [[ $3 == "plus" ]]
then
#  bed="${name}_plus.bed"
  bigwigs=`ls ../../bigwigs_norm/*_for.bw`
  strand="plus"
else
  echo "bad strand"
  exit 1
fi

BS=$2

computeMatrix reference-point -S $bigwigs -R $bed --referencePoint TSS \
	-a 75000 -b 1000 --binSize $BS -p 8 -out ${name}_${strand}.tab.gz \
	--missingDataAsZero --skipZeros 


#optional ags are maxthresh, minthresh, skipzeros, after, before

#for plotprofile, options params are --perGroup, --plotType=fill, --ymin, --ymax, -y, --plotTitle, 
#--plotType=(std,heatmap), 