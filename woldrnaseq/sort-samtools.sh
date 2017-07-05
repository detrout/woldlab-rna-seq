#!/bin/bash

export IN=$1
export OUT=$2
export MEMORY=$3
export THREADS=$4
export PAIRED=$5

if [ -z ${PAIRED} ]; then
    echo "sort single end to ${OUT} (Temp ${_CONDOR_SCRATCH_DIR})"
  cat <( samtools view -H ${IN} ) \
      <( samtools view ${IN}  | \
         sort -S ${MEMORY} -T ${_CONDOR_SCRATCH_DIR} ) | \
      samtools view -@ ${THREADS} -bS - -o ${OUT}
else
    echo "sort paired end to ${OUT} (Temp ${_CONDOR_SCRATCH_DIR})"
  cat <( samtools view -H ${IN} ) \
      <( samtools view ${IN}  | \
         paste -d "|" - - | \
         sort -S ${MEMORY} -T ${_CONDOR_SCRATCH_DIR} | tr '|' '\n' ) | \
      samtools view -@ ${THREADS} -bS - -o ${OUT}
fi
