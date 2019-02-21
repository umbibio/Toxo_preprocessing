#!/bin/bash
# inputs in order 
# 1. Path To Directory of trimed fastq files containing R1s and R2s
# 2. Path To Output Directory to store sam files outputed from this script
# 3. Path To Index file (indexing based name)
# 4. Path To Annotation file
# 5. Path To splice sites file


#IDX="/project/umb_kourosh_zarringhalam/toxoplasmagondii/HS2/Genome/GT1_tran"                                                                                                                                                                
#ANNOTATION="/project/umb_kourosh_zarringhalam/toxoplasmagondii/HS2/Genome/ToxoDB-41_TgondiiGT1.gtf"                                                                                                                                         
#SSIF="/project/umb_kourosh_zarringhalam/toxoplasmagondii/HS2/Genome/ToxoDB-41_TgondiiGT1.ss"       

module load hisat2/2.0.5

TRIM_DIR=$1
SAM_DIR=$2
IDX=$3
ANNOTATION=$4
SSIF=$5


BSUBPARAMS="-n 20 -R rusage[mem=500] -W 02:00 -q short -R span[hosts=1]"
HISATPARAMS="-p 20 -x $IDX --no-mixed --no-discordant -k 1 -X 500 --score-min L,0,-0.6"

for f in $(ls $TRIM_DIR | grep "_val_1" | sort); do
        
    NAME=$(echo $f | cut -f 1 -d '.' | sed 's/_val_1//g' | sed 's/_R1//g')
    JOBID=$NAME".jid"
    READ1=$f
    READ2=$(echo $f | sed 's/_val_1/_val_2/g' | sed 's/_R1/_R2/g')
    #echo $READ1
    #echo $READ2
    CMD="bsub $BSUBPARAMS -J $JOBID -oo hisat2_stdout/$NAME.txt hisat2 $HISATPARAMS -1 $READ1 -2 $READ2 -S $SAM_DIR$NAME.sam"
    echo $CMD
    echo $CMD >> cmd_history.txt

done
	