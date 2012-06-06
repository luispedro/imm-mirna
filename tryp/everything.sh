#!/bin/sh

# Trypanosoma Brucei Genome (Required)
# Chromosomes from http://www.genedb.org/Homepage/Tbruceibrucei927
REF_SOURCE1="ref/Tb927_01_v4.embl"
REF_SOURCE2="ref/Tb927_02_v4.embl"
REF_SOURCE3="ref/Tb927_03_v4.embl"
REF_SOURCE4="ref/Tb927_04_v4.embl"
REF_SOURCE5="ref/Tb927_05_v5.embl"
REF_SOURCE6="ref/Tb927_06_v4.embl"
REF_SOURCE7="ref/Tb927_07_v4.embl"
REF_SOURCE8="ref/Tb927_08_v4.embl"
REF_SOURCE9="ref/Tb927_09_v5.embl"
REF_SOURCE10="ref/Tb927_10_v5.embl"
REF_SOURCE11="ref/Tb927_11_v5.embl"

# TAPyR executable (Required)
CMD_TAPYR="/home/daniel/tapyr13beta/tapyr"

# Reads files (Required)
READS1_TXT="50HA_tab_dl_prefilter.txt"
READS2_TXT="51TY1_tab_dl_prefilter.txt"

REF_BASENAME="ref/Tb927"

# Create fasta from embl chromossomes.
#
# Generates: 
# 	$REF_BASENAME.fa 
create_fasta_from_embl() {
	#echo ">Tb927 Trypanosoma Brucei"
	for i in 1 2 3 4 5 6 7 8 9 10 11 ; do
		echo ">Tb927 Trypanosoma Brucei (Chromosome "$i")"
		cat $(eval echo "\${REF_SOURCE$i}") | egrep "^  " | awk '{ for (i=1; i< NF; ++i) printf("%s", $i); printf("\n") }' 
	done
}
#create_fasta_from_embl > $REF_BASENAME.fa

# Create TAPyR index 
#
# Generates:
#       $REF_BASENAME.fmi
#$CMD_TAPYR I $REF_BASENAME.fa

# Convert reads files to fasta
#
# Generates: 
#	$READS1_TXT.fa
#	$READS2_TXT.fa
convert_reads_txt_to_fasta() {
	cat $READS1_TXT | awk '{ printf(">%s\n%s\n", $1, $2) }' > $READS1_TXT.fa
	cat $READS2_TXT | awk '{ printf(">%s\n%s\n", $1, $2) }' > $READS2_TXT.fa
}
#convert_reads_txt_to_fasta

# Map reads
#
# Generates:
#	$READS1_TXT.sam
#$CMD_TAPYR $REF_BASENAME.fmi $READS1_TXT.fa -i 94
#$CMD_TAPYR $REF_BASENAME.fmi $READS2_TXT.fa -i 94

# Generate counts and concensus
#$CMD_TAPYR c $REF_BASENAME.fa $READS1_TXT.sam

