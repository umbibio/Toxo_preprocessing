#!/bin/bash

module load hisat2/2.0.5

IN_SS=$1
IN_EX=$2
IN_FA=$3
OUT_DIR=$4

bsub -n 6 -R "rusage[mem=4000]" -R "span[ptile=12]" -W 3:59 -q short -J build_index -oo hisat2_stdout/build_index.txt hisat2-build -p 6 --ss $IN_SS --exon $IN_EX $IN_FA $OUT_DIR