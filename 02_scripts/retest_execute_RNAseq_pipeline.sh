#!/usr/bin/env bash
#SBATCH --job-name=TEST_SCRATCH
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=1:00:00
#SBATCH --output=log_test_scratch_RNAseq_pipe_%j.txt


module purge
ml singularity

metadatafile=$1

##execute the RNA-seq_pipeline
#bash summitRNAseq_PE_analyzer.sh /pl/active/onishimura_lab/PROJECTS/RR_ARPE_DELUCA_COLLAB/01_input/RR_ARPE_Deluca_Collab_manifest.txt $SLURM_NTASKS
bash scratch_summitRNAseq_PE_analyzer.sh $metadatafile $SLURM_NTASKS
   ######### MODIFY the SECOND argument to point to YOUR metadata.file ######### 

## OR, you can use a python script
#python RNAseq_analyzer_181011.py ../01_input/metadata_aceticAcid_subset.txt $SLURM_NTASKS


## clean up by zipping .fastq files and deleting extra files
#bash RNAseq_cleanup_mouse_181011.sh /pl/active/onishimura_lab/PROJECTS/RR_ARPE_DELUCA_COLLAB/01_input/RR_ARPE_Deluca_Collab_manifest.txt
   ######### modify the SECOND argument to point to YOUR metadata.file ######### 
