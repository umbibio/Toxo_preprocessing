#!/bin/bash

module load samtools/1.4.1


SAM_DIR=$1

BSUBPARAMS="-n 20 -R rusage[mem=800] -W 00:40 -q short -R span[hosts=1]"

for f in $(ls $SAM_DIR | grep "sam" | grep "RH" | uniq); do
    NAME=$(echo $f | cut -f 1 -d '.')
    JOBID=$NAME".sam.jid"
    #rm -f $SAM_DIR$NAME".bam"
    CMD="bsub $BSUBPARAMS -J $JOBID -oo samview_stdout/$NAME.txt samtools view -@ 20 -o $SAM_DIR$NAME.bam -b $SAM_DIR$NAME.sam"
    $CMD
    #echo $CMD
    echo $CMD >> cmd_history.txt
done
	