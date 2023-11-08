#!/bin/sh --login
#SBATCH -J vcf_table
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64g
#SBATCH --time=24:00:00
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/combined_gvcf/variant_to_table_%j

cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/combined_gvcf/variant_filter
ml purge
ml GCCcore/10.2.0
ml GATK/4.2.0.0-Java-11
ml picard/2.25.0-Java-11

$EBROOTGATK/gatk --java-options "-Xmx60g" VariantsToTable \
-V /mnt/home/goeckeri/SnpEff/snpEff/ALL_filtered_variants_4A_EFF.vcf \
-O ALL_filtered_variants_4A_EFF_table \
-F CHROM \
-F POS \
-F TYPE \
-F REF \
-F ALT \
-F FILTER \
-F QUAL \
-GF GT \
-GF AD \
-GF DP \
-GF PL \
-GF GQ \
-F EFF \
-F ANN \
-F AC \
-F AF

#you need to add the EFF field by creating a database with snpEff and annotating the effects of the variants. 

