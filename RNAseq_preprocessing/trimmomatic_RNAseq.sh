#!/bin/sh --login
#SBATCH -J TRIM_RNAseq
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8g
#SBATCH --time=24:00:00
#SBATCH -o /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/trimmomatic/trimmomatic1_%j
#SBATCH -a 1-48
#SBATCH --export=INFILE=/mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/files_R1_list.txt

module purge
module load Trimmomatic/0.39-Java-11

cd /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/trimmomatic

FILE=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${FILE}
Read1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${Read1}
Read2=$(echo ${Read1}| sed "s/L001_R1_001/L001_R2_001/")
echo ${Read2}
Name=$(echo ${FILE}| sed 's/\/mnt\/scratch\/goeckeri\/RNAseq_bloom\/20230210_mRNASeq_PE150\///' | sed "s/L001_R1_.*//")
echo ${Name}

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE -phred33 -threads 16 \
${Read1} \
${Read2} \
${Name}_R1_paired_trimmed.fastq.gz ${Name}_R1_unpaired_trimmed.fastq.gz \
${Name}_R2_paired_trimmed.fastq.gz ${Name}_R2_unpaired_trimmed.fastq.gz \
ILLUMINACLIP:/mnt/home/goeckeri/genome_files/mont/RNA_seq/TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 \
LEADING:5 \
TRAILING:5 \
MINLEN:45 \
AVGQUAL:15