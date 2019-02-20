#!/bin/bash

module load hisat2/2.0.5
module load python/2.7.14

IN_GTF=$1
OUT_SS=$2

bsub -n 20 -R "rusage[mem=2048]" -W 0:50 -q short -R span[hosts=1] -J splice -oo hisat2_stdout/splice_sites.txt "extract_splice_sites.py $IN_GTF > $OUT_SS" 2>/dev/null