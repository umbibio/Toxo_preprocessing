#!/bin/bash

module load samtools/1.4.1


SAM_DIR=$1

BSUBPARAMS="-n 4 -R rusage[mem=16000] -W 00:50 -q short -R span[hosts=1]"

for f in $(ls $SAM_DIR | grep "bam" | grep "RH" | uniq); do
    NAME=$(echo $f | cut -f 1 -d '.')
    JOBID=$NAME".sort.jid"
    #rm -f $SAM_DIR$NAME".sorted.bam"
    CMD="bsub $BSUBPARAMS -J $JOBID -oo samsort_stdout/$NAME.txt samtools sort -@ 20 -m 1G -o $SAM_DIR$NAME.sorted.bam $SAM_DIR$NAME.bam"
    $CMD
    #echo $CMD
    echo $CMD >> cmd_history.txt
done
