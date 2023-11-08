#!/bin/sh --login
#SBATCH -J stringtie1_subgenomeA_relaxed_gene_level
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=12g
#SBATCH --time=24:00:00
#SBATCH -o /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/stringtie/transcriptome_subgenomeA_relaxed_gene_level_counts_round2QTL_%j
#SBATCH -a 1-48
#SBATCH --export=INFILE1=/mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/STAR/subgenomeA/subgnomeA_relaxed/subgenomeA_relaxed_bam.list

ml purge
ml GCC/8.3.0
ml StringTie/2.1.3

FILE1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE1}`
echo ${FILE1}
Name=$(echo ${FILE1}| sed 's/\/mnt\/scratch\/goeckeri\/RNAseq_bloom\/20230210_mRNASeq_PE150\/STAR\/subgenomeA\/subgenomeA_relaxed\///' | sed "s/_S[0-9][0-9]_subA_relaxed_sorted.bam//")
echo ${Name}

#requantifying after once again fixing a few genes in the QTL region. 

gff_subgenomeA=/mnt/scratch/goeckeri/subgenomeA_updated_June30.gff3

cd /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/stringtie/gene_level_subA_QTL_update_round2/

stringtie -e -p 6 -G $gff_subgenomeA -o ${Name}_QTL_round2.gtf ${FILE1}

#what the fuck is wrong with these output directories. Fuckin hell. I could NOT seem to get these to output in a normal directory. So I just let them show up in the same directory where the file list was and then moved all the .gtfs


#step 1 - make a gtf for alllllll your libraries with the bams.
#step 2 - merge all gtfs for all libraries into one big-ol file that has all non-redundant transcripts in it
#step 3 - run stringtie again with the merged gtf (-G) in expression estimation mode (-e).

#edit: those three steps would be required if i wanted to do isoform-level DE. But after considering the amount of work that would entail
#if we wanted to combine expression of syntelogs for each genotype, I'm not so sure that's where I want to go.

#that being said, I decided to look at ~100 or so genes in the 9 - 10ish Mb region on chr4A to see if doing isoform-level expression would really give us significantly more information than 
#what a simpler gene-level count would. There certainly were some transcripts assembled that didn't seem to be associated with a gene model, maximally 16 - but most of those were single exon or
#very small. I'd say 4-5 actually were convincing of a missing gene annotation. So I'm opting to go the route of gene-level counts. I'll settle with 90%+ gene models being annotated. 


#of course, there's always the chance that a certain isoform of a an annotated gene is expressed in earlies vs lates. But I think I'm not going to worry about that for now...
#I plan to examine the genes expressed in this region pretty throughly anyway. 



#http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

#we'll be comparing syntelogs eventually
#default penalties for multimappers. 
# -u :Turn off multi-mapping correction. In the default case this correction is enabled, and each read that is mapped in n places only contributes 1/n to the transcript coverage instead of 1


#-A flag isn't needed, we never use that table file. counts are derived from the gtf using prepDE.py