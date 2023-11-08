#!/bin/sh --login
#SBATCH -J TRIM_DNAseq
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8g
#SBATCH --time=9:00:00
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/trimmomatic/trimmomatic1_%j
#SBATCH -a 1-14
#SBATCH --export=INFILE=/mnt/home/goeckeri/tart_cherry_bloom/20180803_DNASeq_PE150/DNAseq_R1_only.txt

module purge
module load Trimmomatic/0.39-Java-11
#note - we are not including montmorency in this pipeline because while we know it doesn't have the k haplotype in this form, we don't know whether or not it carries some variant of 'the gene' affecting bloom time. 
#BAL = Balaton = mom
#SF = Surefire = dad

cd /mnt/scratch/goeckeri/pop4_variant_calling/trimmomatic

FILE=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${FILE}
Read1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${Read1}
Read2=$(echo ${Read1}| sed "s/R1/R2/")
echo ${Read2}
Name=$(echo ${FILE}| sed 's/\/mnt\/home\/goeckeri\/tart_cherry_bloom\/20180803_DNASeq_PE150\///' | sed "s/_R1_.*//")
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


