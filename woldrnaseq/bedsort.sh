#!/bin/sh
set -e

export OUTDIR=$1

./bedSort Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str1.out.bg
./bedSort Signal.Unique.str1.out.bg Signal.Unique.str1.out.bg
