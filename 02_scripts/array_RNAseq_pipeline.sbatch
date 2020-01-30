#!/usr/bin/env bash
#SBATCH --job-name=array_rnaseq
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=4:30:00
#SBATCH --output=%A_%a_array_rnaseq.txt
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
echo "\$1 is $1 is $METAFILE"
metadatafile=$(mktemp)
echo "created $metadatafile"
cmd="head -n $SLURM_ARRAY_TASK_ID $METAFILE | tail -1 > $metadatafile"
echo $cmd
eval $cmd

# analyze script
#bash summitRNAseq_PE_analyzer.sh $metadatafile $SLURM_NTASKS
bash safehisat_summitRNAseq_PE_analyzer.sh $metadatafile $SLURM_NTASKS
job_status=$?
if [ $job_status -ne 0 ]
then

    echo "$0: JOB FAILED"
    exit $job_status
    
fi

# run script
#bash RNAseq_cleanup_191204.sh $metadatafile 

rm $metadatafile
