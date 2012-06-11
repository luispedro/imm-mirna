#!/bin/sh

##############################################################

# TAPyR executable
CMD_TAPYR="../tapyr13beta/tapyr"

# FastQC executable
CMD_FASTQC="../FastQC/fastqc"

# DynamicTrim executable (http://solexaqa.sourceforge.net/)
CMD_DYNAMICTRIM="../solexaqa/DynamicTrim.pl"

# LengthSort executable (http://solexaqa.sourceforge.net/)
CMD_LENGTHSORT="../solexaqa/LengthSort.pl"

# DESeq R Script (Requires 'R' to be installed) 
DESEQ_R="mus-deseq.R"

# Mus Musculus genome reference
REF="ref/Mus_musculus.NCBIM37.67.ncrna.fa"

# Reads files (FQ)
READS1="Control-3T3.fq"
READS2="Stimulated_3T3.fq"

##############################################################

REF_FILTERED_BASENAME="Genome.filtered"

READS1_FILTERED=$READS1.filtered.fq
READS2_FILTERED=$READS2.filtered.fq

READS1_TRIMMED=$READS1_FILTERED.trimmed
READS2_TRIMMED=$READS2_FILTERED.trimmed

READS1_SORTED=$READS1_TRIMMED.single
READS2_SORTED=$READS2_TRIMMED.single

READS1_SAM=$READS1_SORTED.sam
READS2_SAM=$READS2_SORTED.sam

READS1_COUNTS=$READS1_SORTED.counts.txt
READS2_COUNTS=$READS2_SORTED.counts.txt

DE_OUTPUT="mus-de-binom.txt"

__preprocess_read() {
	rfile=$1
	rfile_filtered=$2

	echo "Generating FastQC report for "$rfile
	$CMD_FASTQC $rfile --noextract
	
	# get list of overrepresented sequences
	overrep=$(unzip -p ${rfile}_fastqc.zip ${rfile}_fastqc/fastqc_data.txt | sed -n "/Overrepresented/,/>>END_MODULE/p" | head -n-1 | tail -n+3 | grep -v "No Hit")
	
	echo "Removing known overrepresented sequences from "$rfile":"
	echo "$(echo "$overrep" | cut -f4 | sed 's/^/	/')"

	seqs=$(echo "$overrep" | cut -f1 | tr '\n' '|' | sed s/.$//)
	escseqs=$(echo $seqs | sed 's/|/\\|/g')

	sed -n "h;n; /$escseqs/!{x;p;x;p;n;p;n;p;}; /$escseqs/{n;n;}" $rfile > $rfile_filtered

	# tests
	echo "Testing..."
	target=$(echo "$overrep" | awk '{ sum += $2; } END { print sum; }')
	n1=$(wc -l $rfile | cut -f1 -d' ')
	n2=$(wc -l $rfile_filtered | cut -f1 -d' ')
	echo "$( expr $(expr $n1 - $n2) / 4 ) of $target unwanted entries removed."
	echo $(egrep "$seqs" $rfile_filtered | wc -l)" unwanted entries still present."
}

##############################################################

# Preprocess reads files by running FastQC on them and removing
# known overrepresented sequences
#	-> $READS1_FILTERED
#	-> $READS2_FILTERED
preprocess_reads() {
	__preprocess_read $READS1 $READS1_FILTERED
	__preprocess_read $READS2 $READS2_FILTERED
}

# Trim sequences in reads files based on phred score. 
# Sequences will be trimmed to the longest continuous segment
# for which all bases have at minimum the phred score provided.
#	-> $READS1_TRIMMED
#	-> $READS2_TRIMMED
trim_reads() {
	phred=$1 	# Minimum phred score

	echo "Trimming sequences. Phred = "$phred
	$CMD_DYNAMICTRIM $READS1_FILTERED -h $phred
	$CMD_DYNAMICTRIM $READS2_FILTERED -h $phred
}

# Discard sequences smaller than the minimum length provided.
#	-> $READS1_SORTED
#	-> $READS2_SORTED
sort_discard_reads() {
	length=$1	# Minimum sequence length

	echo "Sorting sequences and discarding < "$length"bp"
	$CMD_LENGTHSORT $READS1_TRIMMED -l $length
	$CMD_LENGTHSORT $READS2_TRIMMED -l $length
}

# Cleanup reference genome by removing non-microRNA sequences
#	-> $REF_FILTERED_BASENAME.fa
cleanup_reference() {
	echo "Removing non-microRNA entries from reference."
	cat $REF | sed -n ':i /miRNA/{p; :a n; /^>/!{p;ba}; /^>/{bi} }' > $REF_FILTERED_BASENAME.fa 
}

# Generate TAPyR index 
#	-> $REF_FILTERED_BASENAME.fmi
generate_index() {
	echo "Generating TAPyR index."
	$CMD_TAPYR I $REF_FILTERED_BASENAME.fa
}

# Align reads
#	-> $READS1_SAM
#	-> $READS2_SAM
align_reads() {
	precision=$1	# Requested alignment precision (100 = perfect alignement)

	echo "Aligning reads..."
	$CMD_TAPYR $REF_FILTERED_BASENAME.fmi $READS1_SORTED -i $precision -o $READS1_SAM
	$CMD_TAPYR $REF_FILTERED_BASENAME.fmi $READS2_SORTED -i $precision -o $READS2_SAM
}

# Generate counts
generate_counts() {
	sam=$1 		# SAM file to be processed
	out=$2		# output file

	echo "Generating counts for "$sam
	ns=$(cat $REF_FILTERED_BASENAME.fa | grep "^>" -c)
	cat $sam | grep -v "^@SQ" | cut -f3 | sed -n 's/^R//p' | awk "{ ++cts[\$1] } END { for (i=0; i<$ns; ++i) { if (cts[i] == \"\") cts[i] = 0; print i, cts[i] } }" > $out

	# cat Control-3T3.fq.filtered.fq.trimmed.single.sam | grep "^@SQ" | sed -e 's/^.*SP:\(.*\)\tUR.*$/\1/'
}

# Perform Differential Expression tests.
perform_de() {
	echo "Perform Differential Expression tests."

	# call R script 
	# 	-> mus-de-binom-table.txt
	R CMD BATCH $DESEQ_R

	cat mus-de-binom-table.txt | cut -f8 | tail -n+2 | paste Control-3T3.fq.filtered.fq.trimmed.single.counts.txt Stimulated_3T3.fq.filtered.fq.trimmed.single.counts.txt - | awk 'BEGIN {print "ID      Ctrl    Stim    pval"} { printf("%-8s%-8s%-8s%-8s\n", $1, $2, $4, $5) }' > "mus-de-binom-compare.txt"
}

##############################################################

#preprocess_reads
#trim_reads 30
#sort_discard_reads 16

#cleanup_reference 
#generate_index

#align_reads 100

generate_counts $READS1_SAM $READS1_COUNTS
generate_counts $READS2_SAM $READS2_COUNTS

perform_de
