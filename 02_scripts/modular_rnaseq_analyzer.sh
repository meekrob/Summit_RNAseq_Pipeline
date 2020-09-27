#!/usr/bin/env bash
#SBATCH --job-name=mod_rnaseq
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=1:00:00
#SBATCH --output=%j_mod_rnaseq.txt

# modular_rnaseq_analyzer.sh
# modular_rnaseq_analyzer.sh metadatafile jobstep[s]
filename=$1
shift

# These are the individual functions in the pipeline, defined below
jobsteps=( qc align reformat count )
# This is the step or steps specified on the command line, after the metadata file
jobnames=( $@ )

function parse_and_launch() 
# script how to parse the file and 
# formulate the input, output files based on each jobstep and
# the metadata file
{
    linecount=1 # 1-based
    while read line
    do
        echo "parse_and_launch"
        echo "line: $linecount >> $line <<"
        # ARRAY MODE: only process the $SLURM_ARRAY_TASK_ID line of the arg file
        # skip the others
        if [ -n "$SLURM_ARRAY_TASK_ID" ] && [ $((linecount++)) -ne $SLURM_ARRAY_TASK_ID ]
        then
            continue
        fi

        # Parse metadatafile
        # column 1: left_reads_fastq
        # column 2: right_reads_fastq
        # column 3: name
        read -r left_reads_fastq right_reads_fastq name <<< "$line"

        # define file variables for every step. 
        # They should only differ between input files, not between jobsteps
        
        # These are globally available

        run_fxns # run the functions defined below according to run mode

        ((linecount++))
    
    done < $filename
}

# Define a function for each jobstep as "do_jobstep() {}"
function do_qc()
# PROGRAMS USED:
# FILES USED:
# FILES GENERATED:
{
    echo "in do_qc"
}

slurm_align_args="--ntasks=1 --time=1:00:00"

function do_align()
# PROGRAMS USED:
# FILES USED:
# FILES GENERATED:
{
    echo "in do_align"
}

function do_reformat()
# PROGRAMS USED:
# FILES USED:
# FILES GENERATED:
{
    echo "in do_reformat"
}

function do_count()
# FILES USED:
# PROGRAMS USED:
# FILES GENERATED:
{
    echo "in do_count"
}

# determine run mode and execute desired functions
#

function run_fxns()
{
    
    dependancy_id=""
    for jobstep in ${jobnames[*]}
    do
        fname="do_$jobstep"
        slurm_varname="\$slurm_${jobstep}_args"
        slurm_vals=$(eval echo $slurm_varname)
        if [ -n "$dependancy_id" ]
        then
            slurm_vals="$slurm_vals $dependancy_id"
        fi

        if hash $fname 2>/dev/null
        then
            echo $slurm_vals
            if [ -z "$SLURM_TASK_ID" ] # YOU ARE NOT INSIDE SLURM, LAUNCH WITH SBATCH
            then
                dependancy_id=$(sbatch -b $slurm_vals $THIS_SCRIPT $jobstep)
            else                       # YOU ARE INSIDE SLURM, RUN THE FUNCTION
                cmd="$fname $index"
                echo $cmd
                eval $cmd
            fi
        else
            echo "can't execute $fname, it is not defined"
        fi
    done
}

### Now everything is defined and global variables will determine the requested functions
parse_and_launch
###
