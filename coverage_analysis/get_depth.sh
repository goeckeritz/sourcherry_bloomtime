#!/bin/sh --login
#SBATCH -J bwa_mem_array
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8g
#SBATCH --time=9:59:59
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/coverage_%j
#SBATCH -a 1-14
#SBATCH --export=INFILE=/mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/sam_files.txt

ml GCCcore/9.3.0
ml SAMtools/1.11

cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/coverage/

FILE=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${FILE}
sam=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${sam}
Name=$(echo ${sam}| sed "s/\/mnt\/scratch\/goeckeri\/pop4_variant_calling\/bwa_mem\///" | sed "s/.sam//" | sed "s/_RGs//")
echo ${Name}

#samtools sort ../${Name}_removed_dup.bam > ../${Name}_sorted_dups_removed.bam #DO THIS, MAKE SURE THE SAMPLE BAMS ARE SORTED!!
samtools depth -aa -Q 20 -q 20 /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/${Name}_sorted_dups_removed.bam > ${Name}.depth.post.dups

#these files can be HUUUUUUUGE. And I worked with a small genome... godspeed, friend. 
#The next script you'll want to use is the coverage_loops.sh to do some data reduction.
#Then it's off to making pretty plots. 