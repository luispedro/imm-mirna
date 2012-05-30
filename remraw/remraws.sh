#!/bin/bash

primfile="prims.txt"

seqs=$(cat $primfile | cut -f1 | tr '\n' ' ')
nums=$(cat $primfile | cut -f2 | tr '\n' '+' | sed 's/+/ + /g')

target=$(expr ${nums: 0 : -3})

if test $# -lt 2 ; then
	echo "Syntax: remraws.sh <input filename> <output filename>"
	exit
fi

infile=$1
outfile=$2

cat $infile | python remraw.py $seqs > $outfile 

# check number of entries removed
n1=$(wc -l $infile | cut -f1 -d' ')
n2=$(wc -l $outfile | cut -f1 -d' ')
echo "$( expr $(expr $n1 - $n2) / 4 ) of $target entries removed."

# check if primers are still present in output file
seqspat=$(echo $seqs | tr ' ' '|')
echo "$(egrep $seqspat $outfile -c) unwanted entries found."
