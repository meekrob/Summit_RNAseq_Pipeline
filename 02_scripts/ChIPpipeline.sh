#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=0:35:00
#SBATCH --partition=shas
# --open-mode=append assures that the log files will be appended to from different jobs
# These directives will serve as defaults when submitted via sbatch
# but are comments when run via bash
NTHREADS=${SLURM_NTASKS} # passes --ntasks set above
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
SUBMIT=$0
jobsteps="BWA BAM SPP" #IDR
# the BWA indexes 
bwa_genome=/projects/dcking@colostate.edu/support_data/bwa-index/ce11.unmasked.fa
# PROJECT ORGANIZATION
input_dir=01_FASTQ
align_dir=02_ALIGN
spp_dir=03_SPP
idr_dir=04_IDR

for dir in $align_dir $spp_dir $idr_dir
do
    mkdir -pv $dir
done

# Helper functions
errecho()
{
    1>&2 echo $@
}
run()
{
    echo "running $@"
    time eval $@
}
sb()
{
    sbatch --parsable $@
}
deps() 
{
    vars="$@"
    vars=${vars// /,}
    d=""
    [ -n "$vars" ] && d="-d afterok:$vars"
    echo $d
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
        stages_jids=""
        # skip comment lines
        if [ ${label:0:1} == '#' ] 
        then
            errecho "skipping: $label $rep1_fastq $rep2_fastq $input_fastq" # essentially print 'line' back out
            continue
        fi

        BASE_LOGFILE="$label.jobs.log.$(date +%y%m%d.%H%M)"
        date > $BASE_LOGFILE
        echo -e "\n-----------------------------------------------" >> $BASE_LOGFILE
        echo -e "${label}\t${rep1_fastq}\t${rep2_fastq}\t${input_fastq}" >> $BASE_LOGFILE

        #ALIGN
        # FILENAMES
        rep1_sam=${rep1_fastq/.fastq/.sam}
        rep2_sam=${rep2_fastq/.fastq/.sam}
        input_sam=${input_fastq/.fastq/.sam}
        # JOBS
        if [[ " $jobsteps " =~ " BWA " ]] 
        then
            # you can configure specific time, or other resources, here
            T="--time=0:05:00"
            align_jid1=$(sb $SUBMIT BWA ${input_dir}/$rep1_fastq  ${align_dir}/$rep1_sam)
            align_jid2=$(sb $SUBMIT BWA ${input_dir}/$rep2_fastq  ${align_dir}/$rep2_sam)
            align_jid3=$(sb $SUBMIT BWA ${input_dir}/$input_fastq ${align_dir}/$input_sam)
            stage_jids="$stage_jids $align_jid1 $align_jid2 $align_jid3"
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
            D=$(deps $align_jid1)
            bam_jid1=$(sb $O $D $SUBMIT BAM ${align_dir}/$rep1_sam ${align_dir}/$rep1_bam)
            bam_jid2=$(sb $O $D $SUBMIT BAM ${align_dir}/$rep2_sam ${align_dir}/$rep2_bam)
            bam_jid3=$(sb $O $D $SUBMIT BAM ${align_dir}/$input_sam ${align_dir}/$input_bam)
            stage_jids="$stage_jids $bam_jid1 $bam_jid2 $bam_jid3"
        fi

        # COMPUTE SIGNAL FILES
        # filenames
        # jobs
        

        # PEAK CALLS
        # The OUTPUT FILENAMES are not specified directly, but have the following format 
        # (like L1_1_VS_L1_input.regionPeak.gz), 
        rep1_regionPeak=${rep1_bam%%.bam}_VS_${input_bam%%.bam}.regionPeak.gz
        rep2_regionPeak=${rep2_bam%%.bam}_VS_${input_bam%%.bam}.regionPeak.gz

        # JOBS
        if [[ " $jobsteps " =~ " SPP " ]] 
        then
            D=$(deps $bam_jid1 $bam_jid3)
            ntasks="--ntasks=1"
            spp_jid1=$(sb $ntasks $O $D $SUBMIT SPP $label ${align_dir}/$rep1_bam ${align_dir}/$input_bam)

            D=$(deps $bam_jid2 $bam_jid3)
            spp_jid2=$(sb $ntasks $O $D $SUBMIT SPP $label ${align_dir}/$rep2_bam ${align_dir}/$input_bam)

            stage_jids="$stage_jids $spp_jid1 $spp_jid2"
        fi

        # IDR launch
        # FILENAMES
        idr_out="${label}_1.narrowPeak"
        # JOBS
        if [[ " $jobsteps " =~ " IDR " ]] 
        then
            D=$(deps $spp_jid1 $spp_jid2)
            idr_jid=$(sb $O $D $SUBMIT IDR ${spp_dir}/$rep1_regionPeak ${spp_dir}/$rep2_regionPeak ${idr_dir}/$idr_out)
            stage_jids="$stage_jids $idr_jid"
        fi

        # Use IDR as the rejoin point from all of the branches.
        IDR_JOB_IDS="$IDR_JOB_IDS $idr_jid"

        # Add a command to merge the temporary log files to the base log file
        echo "# To merge the individual jobs logs into this file:" >> $BASE_LOGFILE
        lognames=""
        for jid in $stage_jids
        do
            lognames="$lognames slurm-$jid.out"
        done
        if [ -n "$lognames" ]
        then
            echo "cat $lognames >> $BASE_LOGFILE && rm -v $lognames" >> $BASE_LOGFILE
        fi
        
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
    export TMPDIR=$SLURM_SCRATCH
    export TMP=$TMPDIR
    jobstep=$1
    shift
    errecho "$jobstep SLURM_JOB_ID=$SLURM_JOB_ID"

    source /projects/dcking@colostate.edu/paths.bashrc
#BWA ###################################
    if [ $jobstep == "BWA" ]
    then
        2>&1 bwa | grep ^Program
        2>&1 bwa | grep ^Version
        infastq=$1
        outsam=$2
        cmd="bwa mem -t $SLURM_NTASKS $bwa_genome $infastq > $outsam"
        run $cmd

#BAM ###################################
    elif [ $jobstep == "BAM" ]
    then
        samtools --version
        insam=$1
        outbam=$2
        filter="-F 1536"
        quality="-q 30"
        sort_prefix="$TMPDIR/samsort_$SLURM_JOB_ID"
        cmd="samtools view -@ $SLURM_NTASKS -b $filter $quality -S ${insam} | samtools sort -T $sort_prefix -@ $SLURM_NTASKS -o ${outbam}"
        run $cmd
        cmd="samtools index $outbam"
        run $cmd
        cmd="samtools quickcheck $outbam && rm -v $insam"
        run $cmd

#SPP ###################################
    elif [ $jobstep == "SPP" ]
    then
        prefix=$1
        rep=$2
        input=$3
        outdir=${spp_dir}
        FRAGLEN=150
        SPP=spp # loadbx has bin/spp as wrapper to run_spp.R
        cmd="$SPP -c=$rep   \
             -i=$input \
             -npeak=300000  \
             -odir=${outdir}  \
             -speak=${FRAGLEN} \
             -savr -rf  \
             -savp="${outdir}/$prefix.pdf" \
             -out=${outdir}/$prefix.ccscores"
        run $cmd
    
#IDR ###################################
    elif [ $jobstep == "IDR" ]
    then
        rep1=$1
        rep2=$2
        outfile=$3
        IDR_THRESHOLD=0.05
        cmd="idr --samples $rep1 $rep2 --input-file-type narrowPeak --rank signal.value --soft-idr-threshold ${IDR_THRESHOLD} --plot --use-best-multisummit-IDR"
        run $cmd

#BZ ####################################
    elif [ $jobstep == "BW" ]
    then
        errecho "jobstep BW"

#UNION #################################
    elif [ $jobstep == "UNION" ]
    then
        errecho "jobstep UNION"
    
#NOT DEFINED
    else
        errecho "jobstep $jobstep is not defined. Must be one of $jobsteps"

    fi # END

fi
