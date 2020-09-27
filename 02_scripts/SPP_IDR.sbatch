#!/usr/bin/env bash
#
#
# In order to serve as both the controlling script and individual
# job steps, this script detects whether or not it is being called in
# SLURM job.
SUBMIT=$0
jobsteps="BWA BAM SPP IDR" 
errecho()
{
    1>&2 echo $@
}

run()
{
    errecho "running $@"
    #time eval $@
}
NTASKS=1 # need to change args to pass this in
sb()
{
    errecho "sb: $@"
    while [[ ! $jobsteps =~ $1 ]] # for dry-run, clear potention sbatch args
    do
        shift
    done
    #sbatch -b --ntasks=$NTASKS $@
    jid=$(( ${SLURM_JOB_ID:=0} + 1))
    SLURM_JOB_ID=$jid bash $0 $@
    echo $jid
}

if [ -z "$SLURM_JOB_ID" ]
######## THIS PART OF THE SCRIPT RUNS IN A NORMAL BASH SESSION    ########
######## but launches jobsteps using sbatch (also in this script) ######## 
######## The job steps are invoked by specifying input and output ######## 
######## filenames only. The specific arguments to the programs   ########
########  called within the jobstep definitions in the sbatch     ########
######## block of this script.                                    ########
then
    # Metadata file
    metadatafile=$1
    shift

    [ $# -gt 0 ] && jobsteps="$@"

    # the "filename" portion of each block must always run 
    # in order to make it possible for the pipeline to be 
    # called with an internal starting point, as when
    # the pipeline fails at some specific step.
    while read label rep1_fastq rep2_fastq input_fastq
    do
        # skip comment lines
        if [ ${label:0:1} == '#' ] 
        then
            errecho "skipping: $label $rep1_fastq $rep2_fastq $input_fastq" # essentially print 'line' back out
            continue
        fi

        errecho -e "\n-----------------------------------------------"
        errecho "label: $label"
        errecho "-----------------------------------------------"
        errecho "rep1_fastq: $rep1_fastq"
        errecho "rep2_fastq: $rep2_fastq"
        errecho "input_fastq: $input_fastq"
        errecho "Setting up pipeline for $label"

        #ALIGN
        # FILENAMES
        rep1_sam=${rep1_fastq/.fastq/.sam}
        rep2_sam=${rep2_fastq/.fastq/.sam}
        input_sam=${input_fastq/.fastq/.sam}
        # JOBS
        if [[ " $jobsteps " =~ " BWA " ]] 
        then
            # you can configure specific time, or other resources, here
            align_jid1=$(sb $SUBMIT BWA $rep1_fastq  $rep1_sam)
            align_jid2=$(sb $SUBMIT BWA $rep2_fastq  $rep2_sam)
            align_jid3=$(sb $SUBMIT BWA $input_fastq $input_sam)
        fi

        #CONVERT FORMAT
        # FILENAMES
        rep1_bam=${rep1_sam/.sam/.bam}
        rep2_bam=${rep2_sam/.sam/.bam}
        input_bam=${input_sam/.sam/.bam}
        # JOBS
        if [[ " $jobsteps " =~ " BAM " ]] 
        then
            # you can configure specific time, or other resources, here
            bam_jid1=$(sb --dependancy=afterok:$align_jid1 $SUBMIT BAM $rep1_sam $rep1_bam)
            bam_jid2=$(sb --dependancy=afterok:$align_jid2 $SUBMIT BAM $rep2_sam $rep2_bam)
            bam_jid3=$(sb --dependancy=afterok:$align_jid3 $SUBMIT BAM $input_sam $input_bam)
        fi

        # COMPUTE SIGNAL FILES
        # filenames
        # jobs
        

        # PEAK CALLS
        # FILENAMES
        rep1_narrowPeak=${rep1_bam/.bam/.narrowPeak}        
        rep2_narrowPeak=${rep2_bam/.bam/.narrowPeak}        
        # JOBS
        if [[ " $jobsteps " =~ " SPP " ]] 
        then
            # you can configure specific time, or other resources, here
            spp_jid1=$(sb --dependancy=afterok:$bam_jid1,$bam_jid3 $SUBMIT SPP $rep1_bam $input_bam $rep1_narrowPeak)
            spp_jid2=$(sb --dependancy=afterok:$bam_jid2,$bam_jid3 $SUBMIT SPP $rep2_bam $input_bam $rep2_narrowPeak)
        fi

        # IDR
        # FILENAMES
        idr_out="$label.narrowPeak"
        # JOBS
        if [[ " $jobsteps " =~ " IDR " ]] 
        then
            idr_jid=$(sb --dependancy=afterok:$spp_jid1,$spp_jid2 $SUBMIT IDR $rep1_narrowPeak $rep2_narrowPeak $idr_out)
        fi

        # Use IDR as the rejoin point from all of the branches.
        IDR_JOB_IDS="$IDR_JOB_IDS $idr_jid"

    done < $metadatafile

    # MERGE 
    # FILENAMES
    if [[ " $jobsteps " =~ " MERGE " ]] 
    then
        union_jid=""
    fi

else 
######## THIS PART OF THE SCRIPT RUNS INSIDE SLURM, AND IS CALLED   ########
######## FROM THE BASH SESSION.                                     ########
    errecho "SLURM_JOB_ID=$SLURM_JOB_ID"
    jobstep=$1
    shift

    # BWA
    if [ $jobstep == "BWA" ]
    then
        infastq=$1
        outsam=$2
        cmd="bwa mem -t $SLURM_NTASKS -i $infastq -i $outsam"
        run $cmd

    # BAM
    elif [ $jobstep == "BAM" ]
    then
        insam=$1
        outbam=$2
        cmd="samtools view -@ $SLURM_NTASKS -i $insam | samtools sort -@ $SLURM_NTASKS > $outbam"
        run $cmd
        cmd="samtools index $outbam"
        run $cmd

    # SPP
    elif [ $jobstep == "SPP" ]
    then
        rep=$1
        input=$2
        outfile=$3
        cmd="spp $rep $input > $outfile"
        run $cmd
    
    # IDR
    elif [ $jobstep == "IDR" ]
    then
        rep1=$1
        rep2=$2
        cmd="idr $rep1 $rep2"
        run $cmd

    # BZ
    elif [ $jobstep == "BW" ]
    then
        errecho "jobstep BW"
    
    # Not defined
    else
        errecho "jobstep $jobstep is not defined. Must be one of ${jobstep[@]}"

    fi # END

fi
