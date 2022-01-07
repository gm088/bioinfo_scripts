#!/bin/bash

#window about single point beds
if [ -z $2 ]
then
  echo "usage: bedwin.sh fname halfofwinsize"
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

awk -v var="$winsize" '{print $1"\t"$2-var"\t"$3+var"\t"$4"\t"$5"\t"$6}' ${fname} > ${name}_win.bed
