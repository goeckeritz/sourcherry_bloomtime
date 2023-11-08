#!/bin/sh --login
#SBATCH -J subset
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8g
#SBATCH --time=9:59:59
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/subset_sort_index_%j

ml GCCcore/9.3.0
ml SAMtools/1.11

cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/

#extracting our QTL region from the individuals that have a constant dosage of subgenome A 
for i in $(ls -1 *_sorted_dups_removed.bam); do
	samtools view -h $i "Pcer_chr4A:8324991-16386209" | samtools sort -o- > ./subset_4A/$(echo $i | sed 's/_sorted_dups_removed.bam//')_4A.bam
	samtools index ./subset_4A/$(echo $i | sed 's/_sorted_dups_removed.bam//')_4A.bam
done

#then we need to take the two halves of SF and 27-03-01
#1st region
samtools view -h SF_sorted_dups_removed.bam "Pcer_chr4A:8324991-12583308" | samtools sort -o- > ./subset_4A/SF_4A_1.bam
samtools index ./subset_4A/SF_4A_1.bam
samtools view -h 27-03-01_sorted_dups_removed.bam "Pcer_chr4A:8324991-12583308" | samtools sort -o- > ./subset_4A/27-03-01_4A_1.bam
samtools index ./subset_4A/27-03-01_4A_1.bam

#2nd region -- has a bit of overlap with region 1. (100kb)
samtools view -h SF_sorted_dups_removed.bam "Pcer_chr4A:12483308-16386209" | samtools sort -o- > ./subset_4A/SF_4A_2.bam
samtools index ./subset_4A/SF_4A_2.bam
samtools view -h 27-03-01_sorted_dups_removed.bam "Pcer_chr4A:12483308-16386209" | samtools sort -o- > ./subset_4A/27-03-01_4A_2.bam
samtools index ./subset_4A/27-03-01_4A_2.bam

#for structural	variant	calling, we'll pool k subgenome	samples, i subgenome samples, and call surefire	in the last part of the	QTL and	try to pick out	g from i and k.



#Taking a little outside these markers borders to be safe.(500 kb here)
#extracting alignments to sequences on chr4A from 8324991 to 16386209
#ploidy for each sample in this region:
#Balaton: 1x
#Surefire: 2x -> 3x
#27-02-08: 2x
#27-02-19: 3x
#27-02-52: 2x
#27-03-08: 2x
#27-04-12: 3x
#27-02-65: 1x
#27-03-01: 1x -> 2x
#27-03-25: 0x
#27-03-27: 1x
#27-03-28: 1x
#27-03-46: 1x
#27-04-34: 1x
#split is just about at 12533308 :)   -- we'll go 50 kb on either end to create some overlap in the middle and diminish edge effects. (100kb overlap)

