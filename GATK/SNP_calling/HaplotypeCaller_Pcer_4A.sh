#!/bin/sh --login
#SBATCH -J hap_caller
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=36g
#SBATCH --time=9:59:59
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/haplotype_caller_%j
#SBATCH -a 1-15
#SBATCH --export=INFILE=/mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/4A_bams.txt,INFILE2=/mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/ploidy_4A_QTL_region.txt

cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/

genotype=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${genotype}
ploidy=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE2}`
echo ${ploidy}
sample_name=$(echo ${genotype} | sed 's/\/mnt\/scratch\/goeckeri\/pop4_variant_calling\/bwa_mem\/MarkDups\/subset_4A\///')
echo ${sample_name}

ml purge
ml GCCcore/10.2.0
ml GATK/4.2.0.0-Java-11
ml picard/2.25.0-Java-11

$EBROOTGATK/gatk --java-options "-Xmx36g" HaplotypeCaller \
-R /mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/Pcer_montAv5_p99_400k_scaffolded_assembly.fasta \
-I $(echo ${genotype}) \
-ERC GVCF \
-ploidy $(echo ${ploidy}) \
-ped pedigree_file.txt \
-O $(echo ${sample_name}).g.vcf.gz

#be sure you have a .fai file for your reference genome. Do:
#samtools faidx ref.fasta

#and a dictionary file. For the love of god. Found this out AFTER I was dealing with a bunch of errors. Do:

#java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
#-R <reference fasta>

#the ploidy_4A_QTL_region.txt is a simple file with the following format:

#	~/path_to_directory_where_your_subsetted_bams_are/sample1	2
#	~/path_to_directory_where_your_subsetted_bams_are/sample2	2
#	~/path_to_directory_where_your_subsetted_bams_are/sample3	3

#where the second column is the ploidy of that sample