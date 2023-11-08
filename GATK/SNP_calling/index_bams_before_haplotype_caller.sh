#!/bin/sh

ml purge
ml GCCcore/9.3.0
ml SAMtools/1.11

cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/

for i in $(less 4A_bams.txt); do
	samtools sort $i > $(echo $i | sed 's/.bam//')_sorted.bam
	samtools index $(echo $i)_sorted.bam > $(echo $i).bai
	done
	