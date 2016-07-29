#!/bin/bash

export IN=$1
export OUT=$2
export MEMORY=$3
export THREADS=$4
export PAIRED=$5

if [ -z ${PAIRED} ]; then
    echo "single"
  cat <( samtools view -H ${IN} ) \
      <( samtools view ${IN}  | \
         sort -S ${MEMORY} -T ./ ) | \
      samtools view -@ ${THREADS} -bS - -o ${OUT}
else
    echo "paired"
  cat <( samtools view -H ${IN} ) \
      <( samtools view ${IN}  | \
         paste -d "|" - - | \
         sort -S ${MEMORY} -T ./ | tr '|' '\n' ) | \
      samtools view -@ ${THREADS} -bS - -o ${OUT}
fi
