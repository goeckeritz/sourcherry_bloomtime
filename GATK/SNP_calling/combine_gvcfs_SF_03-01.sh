#!/bin/sh --login
#SBATCH -J combine_those_with_g_haplotypes
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=48g
#SBATCH --time=24:00:00
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/SF_03_01_combine_%j

cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/

ml purge
ml GCCcore/10.2.0
ml GATK/4.2.0.0-Java-11
ml picard/2.25.0-Java-11

$EBROOTGATK/gatk --java-options "-Xmx44g" CombineGVCFs \
-R /mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/Pcer_montAv5_p99_400k_scaffolded_assembly.fasta \
-V 27-03-01_4A_1.bam.g.vcf.gz \
-V 27-03-01_4A_2.bam.g.vcf.gz \
-O 27-03-01_4Ac.bam.g.vcf.gz

#$EBROOTGATK/gatk --java-options "-Xmx36g" CombineGVCFs \
#-R /mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/Pcer_montAv5_p99_400k_scaffolded_assembly.fasta \
#-V SF_4A_1.bam.g.vcf.gz \
#-V SF_4A_2.bam.g.vcf.gz \
#-O SF_4Ab.bam.g.vcf.gz

#man this is taking FOREVER. I ended up making a separate job to get Surefire's combining going too. This is where a better understanding of memory allocation would come in handy.. XD