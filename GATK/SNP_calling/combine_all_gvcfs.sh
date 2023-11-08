#!/bin/sh --login
#SBATCH -J combine_ALL
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64g
#SBATCH --time=72:00:00
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/ALL_combine_GVCFs_%j

cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/combined_gvcf/

ml purge
ml GCCcore/10.2.0
ml GATK/4.2.0.0-Java-11
ml picard/2.25.0-Java-11

$EBROOTGATK/gatk --java-options "-Xmx60g" CombineGVCFs \
-R /mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/Pcer_montAv5_p99_400k_scaffolded_assembly.fasta \
-V ../27-02-08_4A.bam.g.vcf.gz \
-V ../27-02-19_4A.bam.g.vcf.gz \
-V ../27-02-52_4A.bam.g.vcf.gz \
-V ../27-02-65_4A.bam.g.vcf.gz \
-V ../27-03-01_4Ac.bam.g.vcf.gz \
-V ../27-03-08_4A.bam.g.vcf.gz \
-V ../27-03-27_4A.bam.g.vcf.gz \
-V ../27-03-28_4A.bam.g.vcf.gz \
-V ../27-03-46_4A.bam.g.vcf.gz \
-V ../27-04-12_4A.bam.g.vcf.gz \
-V ../27-04-34_4A.bam.g.vcf.gz \
-V ../BAL_4A.bam.g.vcf.gz \
-V ../SF_4Ac.bam.g.vcf.gz \
-O ALL_4A.bam.g.vcf.gz

#what the bloody fuck this took like < 5 minutes while combining 03-01 and SF parts took like 2 DAYS!!!
#Was it because I didn't lower the -Xmx gigs down a little bit from what had been asked for?
#apparently the program needs a little over the top for another function to run