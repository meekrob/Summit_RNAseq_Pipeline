#!/usr/bin/env bash
#SBATCH --job-name=bamcoverage
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
NTHREADS=${SLURM_NTASKS} # passes --ntasks set above
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
date # timestamp
index=${SLURM_ARRAY_TASK_ID:-5}
bamfiles=( $@ )
bamfile=${bamfiles[ $index ]}
bwfile=$(basename $bamfile)
bwfile=${bwfile/.bam/.bw}

module purge
module load singularity/3.3.0
container="/projects/rtpw@colostate.edu/container_images/summitRNAseq.simg"
bamCoverage="singularity exec --bind /scratch/summit/$USER:/tmp $container bamCoverage"

cmd="$bamCoverage \
    -b $bamfile \
    -o $bwfile \
    --outFileFormat bigwig \
    -p $SLURM_NTASKS\
    --normalizeUsing CPM --binSize 1"

echo $cmd
[ -n "$SLURM_NTASKS" ] && time eval $cmd 
status=$?
if [ $status != 0 ]
then
    echo "Command failed with $status"
    exit $status
    
fi
