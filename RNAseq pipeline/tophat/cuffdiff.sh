#!/bin/bash
module load cufflinks/2.2.1

contrasts_dir="/project/umb_kourosh_zarringhalam/toxoplasma_rerun/contrasts/"
cuffquant_combined_dir="/project/umb_kourosh_zarringhalam/toxoplasma_rerun/cuffquant_combined/"
cuffdiff_OUTPUT_DIR="/project/umb_kourosh_zarringhalam/toxoplasma_rerun/cuffdiff/"
GTF_FILE="/project/umb_kourosh_zarringhalam/toxoplasma_rerun/bowtie_index_toxoplasma/ToxoDB35.gtf"
#TOXO_FASTA="/project/umb_kourosh_zarringhalam/toxoplasma_rerun/bowtie_index_toxoplasma/ToxoDB35.fa"
job_logs="/project/umb_kourosh_zarringhalam/toxoplasma_rerun/cuffdiff/logs/"
####### case replicates ####################################
cd ${contrasts_dir}
n=1
while IFS= read -r "case_rep$n"; do
  n=$((n + 1))
done < group11.list.txt

echo $case_rep1
echo $case_rep2
echo $case_rep3

######## control replicates ##################################
cd ${contrasts_dir}
n=1
while IFS= read -r "cntrl_rep$n"; do
  n=$((n + 1))
done < group2.list.txt

echo $cntrl_rep1
echo $cntrl_rep2
echo $cntrl_rep3

######## abundances.cxb files input of cuffdiff ###############

case_rep1_cxb=${cuffquant_combined_dir}${case_rep1}_cuffqunt/abundances.cxb
case_rep2_cxb=${cuffquant_combined_dir}${case_rep2}_cuffqunt/abundances.cxb
case_rep3_cxb=${cuffquant_combined_dir}${case_rep3}_cuffqunt/abundances.cxb
echo 
cntrl_rep1_cxb=${cuffquant_combined_dir}${cntrl_rep1}_cuffqunt/abundances.cxb
cntrl_rep2_cxb=${cuffquant_combined_dir}${cntrl_rep2}_cuffqunt/abundances.cxb
cntrl_rep3_cxb=${cuffquant_combined_dir}${cntrl_rep3}_cuffqunt/abundances.cxb
echo 
experiment_name=${cntrl_rep1%_S*}.over.${case_rep1%_S*}
echo $experiment_name
echo 

######## job submission #######################################
bsub -q long -J ${experiment_name} -e ${job_logs}${experiment_name}".err" -o ${job_logs}${experiment_name}".out" -n 8 -R span[hosts=1] -R rusage[mem=7036] -W 24:00 cuffdiff -p 7 -N -u -L ${case_rep1%_S*},${cntrl_rep1%_S*} -o ${cuffdiff_OUTPUT_DIR}${experiment_name} ${GTF_FILE} ${case_rep1_cxb},${case_rep2_cxb},${case_rep3_cxb} ${cntrl_rep1_cxb},${cntrl_rep2_cxb},${cntrl_rep3_cxb}

#cuffdiff -p 7 -N -u -L B2-P11-intra,B2-P84-intra -o cuffdiff_results_ToxoDB35_default_combined/B2-P84-intra-over-B2-P11-intra /project/umb_kourosh_zarringhalam/toxoplasmagondii/Genome/bowtie_index_toxoplasm/ToxoDB35.gtf cuffquant_results_ToxoDB35_default_combined/B2-P11-intracellular_S4_out/abundances.cxb,cuffquant_results_ToxoDB35_default_combined/B2-P11-intra-bio-rep-2_S5_out/abundances.cxb,cuffquant_results_ToxoDB35_default_combined/B2-P16-intra-bio-rep-3_S1_out/abundances.cxb cuffquant_results_ToxoDB35_default_combined/B2-P84-intracellular_S1_out/abundances.cxb,cuffquant_results_ToxoDB35_default_combined/B2-P86-intra-bio-rep-2_S6_out/abundances.cxb,cuffquant_results_ToxoDB35_default_combined/B2-P86-intra-bio-rep-3_S2_out/abundances.cxb