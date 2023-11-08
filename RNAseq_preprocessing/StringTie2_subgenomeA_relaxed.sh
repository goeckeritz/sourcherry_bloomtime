#!/bin/sh --login
#SBATCH -J stringtie_subgenomeA_relaxed
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=12g
#SBATCH --time=24:00:00
#SBATCH -o /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/stringtie/transcriptome_subgenomeA_relaxed_%j
#SBATCH -a 1-48
#SBATCH --export=INFILE1=/mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/STAR/subgenomeA/subgnomeA_relaxed/filtered_BAMs/filtered_sorted_bams.txt

ml purge
ml GCC/8.3.0
ml ml StringTie/2.1.3

FILE1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE1}`
echo ${FILE1}
Name=$(echo ${FILE1}| sed 's/\/mnt\/scratch\/goeckeri\/RNAseq_bloom\/20230210_mRNASeq_PE150\/STAR\/subgenomeA\/subgenomeA_relaxed\/filtered_BAMs\///' | sed "s/_subA_relaxed_filtered_sorted.bam//")
echo ${Name}

subgenomeA=/mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/chr4_QTL_update/subgenomeA.fasta
gff_subgenomeA=/mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/chr4_QTL_update/subgenomeA_chr4_QTL_updated.gff3

cd /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/stringtie/


stringtie -p 14 -G $gff_subgenomeA -t -c 1.5 -f 0.05 -o ${Name}.gtf


#step 1 - make a gtf for alllllll your libraries with the bams.
#step 2 - merge all gtfs for all libraries into one big-ol file that has all non-redundant transcripts in it
#step 3 - run stringtie again with the merged gtf (-G) in expression estimation mode (-e).

#I believe we use that for DESeq2 :) Then the exciting part begins :D



