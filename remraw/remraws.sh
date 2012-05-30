#!/bin/bash

dir_name=$(dirname ${BASH_SOURCE[0]})

if [ $# -lt 3 ]; then
	echo "Syntax: remraws.sh <input filename> <seqs file> <output filename>"
	exit
fi

infile=$1
seqsfile=$2
outfile=$3

seqs=$(cat $seqsfile | cut -f1 | tr '\n' ' ')
nums=$(cat $seqsfile | cut -f2 | tr '\n' '+' | sed 's/+/ + /g')

target=$(expr ${nums: 0 : -3})

cat $infile | python $dir_name/remraw.py $seqs > $outfile 

# check number of entries removed
n1=$(wc -l $infile | cut -f1 -d' ')
n2=$(wc -l $outfile | cut -f1 -d' ')
echo "$( expr $(expr $n1 - $n2) / 4 ) of $target entries removed."

# check if primers are still present in output file
seqspat=$(echo $seqs | tr ' ' '|')
echo "$(egrep $seqspat $outfile -c) unwanted entries found."
