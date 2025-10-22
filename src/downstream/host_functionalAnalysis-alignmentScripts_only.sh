#!/bin/bash -login
#SBATCH -t 04-00:00:00
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jknorris@ucdavis.edu
#SBATCH -p med2
#SBATCH --output=/home/jknorris/.slurm_output/Btaurus/%x.%j.out
#SBATCH -A ctbrowngrp

set -e  # Exit immediately on error
set -o pipefail
set -u
source ~/.bashrc

# Job information logging (improved)
echo "Job ${SLURM_JOB_NAME}_${SLURM_JOB_ID} running on:"
echo "  NODELIST: ${SLURM_NODELIST}"
echo "  Job details:"
scontrol show job $SLURM_JOBID


# Define nproc explicitly for better control and portability
nproc=$SLURM_CPUS_PER_TASK

source ~/.bashrc    
echo "STARTING"
OUTDIR="/home/jknorris/storage/ctbrownroot/Hugo/BOSTAURUS/"
REFERENCE="${OUTDIR}/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa"
GTF="${OUTDIR}/Bos_taurus.ARS-UCD1.3.113.gtf"
HISAT2_INDEX="${OUTDIR}/HISAT2/Bos_taurus.ARS-UCD1.3"
FASTQC_OUTDIR="${OUTDIR}/HISAT2/fastqc"
ALIGN_OUTDIR="${OUTDIR}/HISAT2/ALIGN"
sample_name=$(get_samplename "$1" | sed 's/_S.*//')
output_prefix="${ALIGN_OUTDIR}/${sample_name}"

build_hisat2_index() {
    if [ ! -s ${HISAT2_INDEX}.1.ht2l ]; then
        mkdir -p ${OUTDIR}/HISAT2/
	load_hisat2
        hisat2-build-l -p $(nproc) ${REFERENCE} ${HISAT2_INDEX}
    fi
}

run_fastqc() {
    local input1=$1
    local input2=$2

    [ ! -d ${FASTQC_OUTDIR} ] && mkdir -p ${FASTQC_OUTDIR}

    load_fastqc
	if [[ ! -s "${FASTQC_OUTDIR}/${sample_name}"*_fastqc.zip && ! -s "${FASTQC_OUTDIR}/${sample_name}"*_fastqc.html ]]; then
        	fastqc ${input1} ${input2} -o ${FASTQC_OUTDIR} --threads $(nproc)
	fi
}

run_hisat2() {
local input1=$1
local input2=$2

[ ! -d ${ALIGN_OUTDIR} ] && mkdir -p ${ALIGN_OUTDIR}

if [ ! -s ${output_prefix}.sam ] 
then
    load_hisat2
    hisat2 -q -t --threads $(nproc) --quiet --summary-file ${output_prefix}.summary -x ${HISAT2_INDEX} \
	-1 ${input1} -2 ${input2} -S ${output_prefix}.sam

    load_samtools
    samtools collate -@ $(nproc) -o ${output_prefix}.sort.sam ${output_prefix}.sam
	mv ${output_prefix}.sort.sam ${output_prefix}.sam
    [ ! -s  ${output_prefix}.flagstat ] && samtools flagstat -@ $(nproc) -O tsv ${output_prefix}.sam > ${output_prefix}.flagstat
fi

}

run_featurecounts() {

if [ ! -s ${output_prefix}.counts ] 
then
    load_featurecounts
    featureCounts -p -t exon -g gene_id -T $(nproc) -a ${GTF} -o ${output_prefix}.counts ${output_prefix}.sam
	load_pigz 
	pigz -p $(nproc) ${1}
fi

}

#############################################################################################################################################

process_sam() {
    echo "Processing SAM: ${sample_name}"

        if [[ ${1} =~ ".gz" ]]; then
		load_pigz
            pigz -d -p $(nproc) ${1}
            pigz -d -p $(nproc) ${1/_R1_001_/_R2_001_}
	R1=${1/.gz/}
	R2=${1/_R1_001_trim.fixed.gz/_R2_001_trim.fixed}
        	echo "${sample_name} - pigz done"
	else 	
	R1=${1}
	R2=${1/_R1_001_/_R2_001_}
        fi

if [[ ! "$R1" == *.fq* || ! "$R2" == *.fq* ]]; then
mv ${R1} ${R1}.fq
R1=${R1}.fq
mv ${R2} ${R2}.fq
R2=${R2}.fq
fi

    run_fastqc ${R1} ${R2}
    run_hisat2 ${R1} ${R2}
    run_featurecounts ${output_prefix}.sam
}

#############################################################################################################################################

main() {
build_hisat2_index
process_sam ${1}
}

main "$@"

printf "Time Used: $(squeue -h -j $SLURM_JOBID -o "%M")\nTime Left: $(squeue -h -j $SLURM_JOBID -o "%L")\nTime Requested: $(squeue -h -j $SLURM_JOBID -o "%l")\nDays Until Halloween: $(fHalloween)"
exit 0
