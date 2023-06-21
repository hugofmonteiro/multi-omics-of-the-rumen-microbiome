#! /bin/bash
#
#SBATCH --mail-user=hmonteiro@ucdavis.edu       # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS - ALL, NONE, BEGIN, END, FAIL, REQUEUE
#SBATCH -J Kraken2db                            # JOB ID
#SBATCH -e Kraken2db.j%j.err                    # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o Kraken2db.j%j.out                    # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH -c 8                                    # NUMBER OF PROCESSORS PER TASK
#SBATCH --mem=200Gb                             # MEMORY POOL TO ALL CORES
#SBATCH --time=07-36:00:00                      # REQUESTED WALL TIME
#SBATCH -p med2                                 # PARTITION TO SUBMIT TO

# initialize conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate snakemake
cd ~/Rumen_Microbiome_Genomics/1_Sequences_Guanhui/kraken2-2.1.2/

# fail on weird errors
set -e
set -x

### YOUR COMMANDS GO HERE ###
#build custom database
#kraken2-build --download-taxonomy --db rumen_db
#kraken2-build --download-library bacteria --db rumen_db
#kraken2-build --download-library archaea --db rumen_db
#kraken2-build --download-library viral --db rumen_db
#kraken2-build --download-library fungi --db rumen_db
#kraken2-build --download-library protozoa --db rumen_db
#kraken2-build --download-library plasmid --db rumen_db
#kraken2-build --download-library UniVec_Core --db rumen_db
#kraken2-build --build --threads 8 --db rumen_db
#adding plant genomes to the library
#to be done after we have all feed ingredients used in the studies so we can download the genome of all these plants
#PHiX
kraken2-build --download-taxonomy --db Cow_and_PhiX_db
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz | gunzip > phiX.fa
kraken2-build --add-to-library phiX.fa --db Cow_and_PhiX_db --no-masking
rm phiX.fa
#adding the just downloaded complete microbial_plant db from rumen_complete_db
# 
#adding cow to the library
kraken2-build --add-to-library bosTaurus.fa --db Cow_and_PhiX_db
#at the end build all downloaded libraries
kraken2-build --build --threads 8 --db Cow_and_PhiX_db
kraken2-build --clean --db Cow_and_PhiX_db
### YOUR COMMANDS GO HERE ###


# Print out values of the current jobs SLURM environment variables
env | grep SLURM

# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
