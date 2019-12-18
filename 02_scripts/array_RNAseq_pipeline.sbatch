#!/usr/bin/env bash
#SBATCH --job-name=array_rnaseq
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=0:30:00
#SBATCH --output=log_RNAseq_pipe_%j_%a.txt
## Run this script like:
##
## sbatch --array 1-J array_RNAseq_pipeline.sbatch metadata.txt
##
## Where J is the number of lines in metadata.txt (no headers allowed in that file).

# set tmpfile location
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR

# make a one-line metadata file with one specific line
METAFILE=$1
metadatafile=$(mktemp)
head -n $SLURM_ARRAY_ID $METAFILE | tail -1 > $metadatafile

# change MODE or comment out to run the clean script
MODE="RUN" 

if [ "$MODE" -eq "RUN" ]
then

    bash RNAseq_analyzer_191204.sh $metadatafile $SLURM_NTASKS

else # CLEAN

    bash RNAseq_cleanup_191204.sh $metadatafile 
fi

rm $metadatafile
