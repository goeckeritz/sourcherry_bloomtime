#!/bin/sh --login
#SBATCH -J index_subgenomeA_relaxed
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=8g
#SBATCH --time=1:00:00
#SBATCH -o /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/STAR/subgenomeA/subgnomeA_relaxed/index_subgenomeA_relaxed_%j
#SBATCH -a 1-48
#SBATCH --export=INFILE1=/mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/STAR/subgenomeA/subgnomeA_relaxed/bams.txt

ml purge
ml GCC/6.4.0-2.28  OpenMPI/2.1.1
ml SAMtools/1.9

FILE1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE1}`
echo ${FILE1}

cd /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/STAR/subgenomeA/subgnomeA_relaxed/

samtools index -@ 6 ${FILE1}

