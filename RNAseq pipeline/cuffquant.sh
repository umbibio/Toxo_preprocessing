#!/bin/bash

#this script quantifies expression of genes in every sample and produces a binary cxb file to be used by cuffdiff

################### loading Modules ##############################################
 

module load cufflinks/2.2.1

#################### Job submission parameter ####################################

PROCESSORS=3
MEMORY="49152"
DURATION="24:00"
QUEUE="long"

################## Input Output directories #####################################

GTF_FILE="/project/umb_kourosh_zarringhalam/toxo_new/Genome/bowtie_index_toxoplasma/ToxoDB35.gtf"
SAMPLE_FILE="/project/umb_kourosh_zarringhalam/toxo_new/Genome/FASTQ_FILES/sample.list.txt"
thout_DIR="/project/umb_kourosh_zarringhalam/toxo_new/Genome/tophat_out_run6/"
cuffquant_OUT_DIR="/project/umb_kourosh_zarringhalam/toxo_new/Genome/cuffquant_out_run6/"

for sample in $(cat ${SAMPLE_FILE}); do

cd  ${thout_DIR}${sample}_thout/
BAMFILE=$(ls *accepted_hits.bam*)
bsub -J "cuff_quant_run6" -q long -W 24:00 -n 8 -R "rusage[mem=49152] span[hosts=1]" cuffquant -p 3 -v ${GTF_FILE} ${BAMFILE} -o ${cuffquant_OUT_DIR}${sample}_cuffqunt 

done
 



