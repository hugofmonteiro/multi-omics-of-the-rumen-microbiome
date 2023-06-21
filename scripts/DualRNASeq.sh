#! /bin/bash
#
#SBATCH --mail-user=hmonteiro@ucdavis.edu   # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                     # NOTIFICATIONS OF SLURM JOB STATUS - ALL, NONE, BEGIN, END, FAIL, REQUEUE
#SBATCH -J DualRNASeq                       # JOB ID
#SBATCH -e Sourmash_classify.j%j.err        # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o Sourmash_classify.j%j.out        # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH -c 16                               # NUMBER OF PROCESSORS PER TASK
#SBATCH --mem=200Gb                         # MEMORY POOL TO ALL CORES
#SBATCH --time=07-36:00:00                  # REQUESTED WALL TIME
#SBATCH -p med2                             # PARTITION TO SUBMIT TO

# fail on weird errors
set -e
set -x

# initialize conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate snakemake

# cd to the folder of desire, if necessary

# paste below the code to be used:

snakemake -p -j 16 -s Snakefile_DualRNASeq.sn --use-conda

# Print out values of the current jobs SLURM environment variables
env | grep SLURM

# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
