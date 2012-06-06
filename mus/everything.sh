#!/bin/sh

CMD_TAPYR="../tapyr13beta/tapyr"

REF_PREFILTERED="Mus_musculus.NCBIM37.67.ncrna.fa"
REF_FILTERED_BASENAME="Mus_musculus.NCBIM37.67.ncrna.filtered"

READS1="Control-3T3_Clean.fq.trimmed.single"
READS2="Stimulated_3T3_Clean.fq.trimmed.single"

cleanup_reference() {
	cat $REF_PREFILTERED | python cleanup.py
}
cleanup_reference > $REF_FILTERED_BASENAME.fa

# Generate TAPyR index -> $REF_FILTERED_BASENAME.fmi
$CMD_TAPYR I $REF_FILTERED_BASENAME.fa

# Align reads 100%
$CMD_TAPYR $REF_FILTERED_BASENAME.fmi $READS1 -i 100
