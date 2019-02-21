#!/bin/bash

module load subread/1.6.2

SAM_DIR=$1

ANNOTATION="../Genome/ToxoDB-41_TgondiiGT1.gtf"

BSUBPARAMS="-n 4 -R rusage[mem=16000] -W 3:59 -q short -R span[hosts=1]"

FCPARAMS="-M -O -p -B -C -T 6"

for f in $(ls $SAM_DIR | grep "sorted" | grep "RH" | uniq); do
    NAME=$(echo $f | cut -f 1 -d '.')
    JOBID=$NAME".count.jid"
    CMD="bsub $BSUBPARAMS -J $JOBID -oo feature_count_stdout/$NAME.txt featureCounts $FCPARAMS -a $ANNOTATION -o $SAM_DIR$NAME.counts.txt $SAM_DIR$NAME.sorted.bam "
    $CMD
    echo $CMD >> cmd_history.txt
    #echo $CMD
done
