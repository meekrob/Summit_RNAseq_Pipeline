#!/usr/bin/env bash

#SBATCH --job-name=RWC21_JM149_Embryo_FANS_RNAseq
#SBATCH --nodes=1
#SBATCH --ntasks=24      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas #modify this to reflect which queue you want to use. Options are 'shas' and 'shas-testing'
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=2:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=stdout/log_RNAseq_pipe_%j.txt


module purge
ml singularity

## execute the RNA-seq_pipeline
bash summitRNAseq_PE_analyzer.sh ../01_input/RWC21_metadata.csv $SLURM_NTASKS
   ######### MODIFY the SECOND argument to point to YOUR metadata.file ######### 

## OR, you can use a python script
#python RNAseq_analyzer_181011.py ../01_input/metadata_aceticAcid_subset.txt $SLURM_NTASKS


## clean up by zipping .fastq files and deleting extra files
#bash RNAseq_cleanup_mouse_181011.sh ../04_testing/metadata_mouse.txt
   ######### modify the SECOND argument to point to YOUR metadata.file ######### 
