NAME:
ChIPpipeline.sh

USAGE:
bash ChIPpipeline.sh - show configuration, jobsteps, look for warnings
bash ChIPpipeline.sh metadatafile.txt [ jobsteps ] - run with a metadata control file
                                                   - without jobsteps, all are run
                                                   - specific steps may be specified

ALTERNATE USAGE (advanced):
sbatch ChIPpipeline.sh JOBSTEP stage infile outfile - this will run in SLURM. The pipeline
                                                    - will invoke this mode of the script.
                                                    - If you wish to run a single job with
                                                    - manual arguments (i.e. filenames other
                                                    - than the chain created by metafile.txt),
                                                    - then you can choose this option.
                                                    - It is a good idea to understand how
                                                    - the step is implemented in the script's
                                                    - lower half (SLURM mode) in order to use
                                                    - this.
DESCRIPTION:
    This script is modified from the original modENCODE ChIP-seq analysis pipeline, which is 
separated out into separate stages for each step. This script merges those files into a single 
script and invokes them using dependency arguments to sbatch in order to coordinate branching 
and merging trajectories in the pipeline.

    Also, the merging of developmental stages and aggregation of signal over them differs
from methods used by modENCODE.

REQUIRED PROGRAMS:
    bash > 4.0, SLURM, R
    bwa 
    spp
    UCSC UserApps
    bedtools

REQUIRED FILES:
    bwa genome indexes
    black list BED file.
    A file with the chromosome sizes.


Last edit:
    Mon Sep 28 12:21:23 MDT 2020
