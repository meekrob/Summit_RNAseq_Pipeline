#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=stdout/fastp-%j.out
#SBATCH --qos=normal
#SBATCH --time=00:10:00

module purge

ml singularity


cmd1="singularity exec ../01_container/summitRNAseq.simg fastp -i ../02_input/EO055-guts1_CGATGT_L005_R1_001.fastq.gz -o ../03_fastp_output/testOut.fastq.qz --json=fastp.json --html=fastp.html"

echo "run fastp"
echo $cmd1 # echo command to .out
eval time $cmd1 # run the command and measure the time it takes

cmd2="mv ./fastp.html ../03_fastp_output"
echo $cmd2
eval $cmd2 

cmd3="mv ./fastp.json ../03_fastp_output"
echo $cmd3
eval $cmd3 
