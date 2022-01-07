#!/bin/bash

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
if [ -z $2 ]
then
  echo "usage: bedwin.sh fname size"
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
winsize=$2

check=`awk '{print $3-$2}' $1 | sort | uniq`

if [ "$check" -ne 1 ]
then
  echo "regions need to be of length 1"
  exit -1
fi


awk -F'\t' -vOFS='\t' -v var="$winsize" '{if($6=="+"){$2=$2-var}{print $0}}' $1 \
| awk -F'\t' -vOFS='\t' -v var="$winsize" '{if($6=="-"){$3=$3+var}{print $0}}' | cut -f 1-6 > ${name}_precede.bed




