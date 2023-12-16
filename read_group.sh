#!/bin/bash
# Run with: bash read_group.sh [ref] [file_r1] [file_r2] [sample] [output]

HEADER=$(zcat $2 | head -n 1)
REF=$1
R1=$2
R2=$3
OUT=$5

FLWCL=$(echo $HEADER | cut -d':' -f3)
LN=$(echo $HEADER | cut -d':' -f4)
BARCODE=$(echo $HEADER | cut -d':' -f10)

ID=$(echo $FLWCL.$LN)
PL="ILLUMINA"
PU=$(echo $FLWCL.$LN.$BARCODE)
LB=$BARCODE
SM=$4

bwa-mem2 mem -t 12 -R $(echo "@RG\tID:${ID}\tPL:${PL}\tPU:${PU}\tLB:${LB}\tSM:${SM}") $REF $R1 $R2 | samtools sort -o $OUT