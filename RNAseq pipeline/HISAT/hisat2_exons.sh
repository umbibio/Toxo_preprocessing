#!/bin/bash

module load hisat2/2.0.5
module load python/2.7.14

IN_GTF=$1
OUT_EX=$2

bsub -n 20 -R "rusage[mem=2048]" -W 0:50 -q short -R span[hosts=1] -J exon -oo hisat2_stdout/exon.txt "extract_exons.py $IN_GTF >$OUT_EX" 2>/dev/null