#!/bin/sh --login
#SBATCH -J bwa_mem_array
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8g
#SBATCH --time=9:59:59
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/mark_dups_%j
#SBATCH -a 1-14
#SBATCH --export=INFILE=/mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/sam_files.txt


cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/

ml purge
ml GCCcore/10.2.0
ml GATK/4.2.0.0-Java-11
ml picard/2.25.0-Java-11

#I can use the sam_files.txt again for file names
FILE=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${FILE}
sam=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${sam}
Name=$(echo ${sam}| sed "s/\/mnt\/scratch\/goeckeri\/pop4_variant_calling\/bwa_mem\///" | sed "s/.sam//" | sed "s/_RGs//")
echo ${Name}

java -jar $EBROOTPICARD/picard.jar SortSam \
-I /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/set_tags/${Name}.fixed.bam \
-O /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/set_tags/${Name}_query_sorted.bam \
-SORT_ORDER queryname

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
-I /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/set_tags/${Name}_query_sorted.bam \
-O ${Name}_marked_dup.bam \
-M ${Name}_marked_dup_metrics.txt \
-DS SUM_OF_BASE_QUALITIES

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
-I /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/set_tags/${Name}_query_sorted.bam \
-O ${Name}_removed_dup.bam \
-M ${Name}_removed_dup_metrics.txt \
-DS SUM_OF_BASE_QUALITIES \
--REMOVE_DUPLICATES true
      



#This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA
#The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method)
#Invoking the TAGGING_POLICY option, you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag
#When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.
#Wouldn't that mean we wanna sort by query?
