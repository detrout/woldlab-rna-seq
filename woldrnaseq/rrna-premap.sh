#!/bin/bash
# This script was based on instructions provided
# by Saurabh Agarwal

GENOME_ROOT=$1
GENOME=$2
THREADS=$3
INPUT=$4
OUTPUT=$5

BOWTIE2_INDICES=${GENOME_ROOT}/rrna/BOWTIE2_INDICES
GENOME_FAI=${GENOME_ROOT}/rrna/${GENOME}_plus_rRNA.fa.fai

/woldlab/castor/proj/programs/bowtie2-2.1.0/bowtie2 \
  -q --threads ${THREADS} --time -x ${BOWTIE2_INDICES}/${GENOME}_premap_BOWTIE2 \
  --score-min L,-0.8,-0.8 -U ${INPUT} | \
samtools view -F 0x4 -F 0x100 -F 0x200 | \
samtools view -but ${GENOME_FAI} | \
cut -f 3 | sort | uniq -c > ${OUTPUT}

#save the premapped bamfile
#samtools sort -@ 2 -m 8G - -o ${OUTPUT}
