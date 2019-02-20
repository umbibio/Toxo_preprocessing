#!/bin/bash

module load gffread/0.9.8

IN_GFF=$1
OUT_GTF=$2

bsub -n 20 -R rusage[mem=2048] -W 0:50 -q short -R span[hosts=1] -J gffTogtf -oo hisat2_stdout/gff_to_gtf.txt gffread -T $IN_GFF -o $OUT_GTF
