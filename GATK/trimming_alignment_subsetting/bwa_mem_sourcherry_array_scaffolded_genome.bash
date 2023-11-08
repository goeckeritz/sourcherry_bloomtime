#!/bin/sh --login
#SBATCH -J bwa_mem_array
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8g
#SBATCH --time=9:59:59
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/bwa_mem_%j
#SBATCH -a 1-14
#SBATCH --export=INFILE=/mnt/scratch/goeckeri/pop4_variant_calling/trimmomatic/paired_reads_R1.txt

module load GCC/6.4.0-2.28 OpenMPI/2.1.1
module load BWA/0.7.17


cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/

FILE=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${FILE}
Read1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${Read1}
Read2=$(echo ${Read1}| sed "s/R1/R2/")
echo ${Read2}
ID=$(echo ${FILE}| sed 's/\/mnt\/scratch\/goeckeri\/pop4_variant_calling\/trimmomatic\///' | sed 's/.*_S/S/' | sed 's/_[R1,R2].*//')
echo ${ID}
ILLUMINA=$(echo "ILLUMINA")
SAMPLE=$(echo ${FILE}| sed 's/\/mnt\/scratch\/goeckeri\/pop4_variant_calling\/trimmomatic\///' | sed "s/_S.*//") 


bwa mem -t 20 -R "@RG\tID:${ID}\tSM:${SAMPLE}\tPL:${ILLUMINA}\tLB:${SAMPLE}_1" /mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/Pcer_montAv5_p99_400k_scaffolded_assembly ${Read1} ${Read2} > ${SAMPLE}_RGs.sam

