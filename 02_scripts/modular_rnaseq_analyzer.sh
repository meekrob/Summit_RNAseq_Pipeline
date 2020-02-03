#!/usr/bin/env bash
# SLURM HEADERS
# SLURM HEADERS
# SLURM HEADERS
# SLURM HEADERS
# SLURM HEADERS
# SLURM HEADERS
# SLURM HEADERS
# set -x

filename="metadata.txt"

jobsteps=( qc align reformat count )
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
    
    for jobstep in ${jobnames[*]}
    do
        fname="do_$jobstep"
        echo "doing jobstep: $fname"
        if hash $fname 2>/dev/null
        then
            $fname $index
        else
            echo "can't execute $fname, it is not defined"
        fi
    done
}

### Now everything is defined and global variables will determine the requested functions
parse_and_launch
###
