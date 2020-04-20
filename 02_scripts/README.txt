These scripts test different approaches to parallelize the RNA-seq pipeline
written by Erin Osborne Nishimura, and modified  Robert Williams at
https://github.com/rtpwilliams/Summit_RNAseq_Pipeline.  Currently, only the
array method has been tested.  You must have a metadata file that has no
header. Say it is named "metadata.txt" Then you run the calling script with the
array argument, like so:

    sbatch --array 1-J array_RNAseq_pipeline.sbatch metadata.txt

... where J is the number of lines in metadata.txt.

The script calls safehisat_summitRNAseq_PE_analyzer.sh, so named because of
inconsistencies with running hisat2 in parallel.  It is limited to running
hisat2 with a single thread, but the entirety of the pipeline is run in parallel, 
one per sample, and the other commands can make use of the value of ntasks.

REQUIRED:
    Edit array_RNAseq_pipeline.sbatch to modify the --ntasks parameter in the
slurm header directives.  Also, the value of --time is the allotted time for
each sample, rather than the collective workload of all the samples in the
metadata file.  

    You must edit safehisat_summitRNAseq_PE_analyzer.sh to use the
pertanent information and command paths for your data.  These are variables
set at the top of the script.

MERGING OUTPUT:
The array pipeline will write separate files for each sample during the featureCounts step. To merge them, do the following:
    module load R
    path_to_repo/02_scripts/merge_counts.RScript sample1 sample2 sample3
This will create a file called "counts.txt" with the featureCounts merged from each sample.

