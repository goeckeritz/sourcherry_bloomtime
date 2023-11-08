#!/bin/sh --login
#SBATCH -J stringtie1_subgenomeA_relaxed_gene_level
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=12g
#SBATCH --time=24:00:00
#SBATCH -o /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/stringtie/transcriptome_subgenomeA_relaxed_gene_level_fpkm_%j
#SBATCH -a 1-47
#SBATCH --export=INFILE1=/mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/STAR/subgenomeA/subgnomeA_relaxed/subgenomeA_relaxed_bam_outlier_removed.list

ml purge
ml GCC/8.3.0
ml StringTie/2.1.3

FILE1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE1}`
echo ${FILE1}
Name=$(echo ${FILE1}| sed 's/\/mnt\/scratch\/goeckeri\/RNAseq_bloom\/20230210_mRNASeq_PE150\/STAR\/subgenomeA\/subgenomeA_relaxed\///' | sed "s/_S[0-9][0-9]_subA_relaxed_sorted.bam//")
echo ${Name}

#requantifying after once again fixing a few genes in the QTL region. 

gff_subgenomeA=/mnt/scratch/goeckeri/round2_update/tmp_files/rename_the_genes/subgenomeA_chr4_QTL_update2.gff3

cd /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/stringtie/gene_level_subA_QTL_update_round2/

stringtie -e -p 8 -G $gff_subgenomeA -A ${Name}_QTL_round2_4fpkm -o ${Name}_QTL_round2_4fpkm.gtf ${FILE1}

#what the fuck is wrong with these output directories. Fuckin hell. I could NOT seem to get these to output in a normal directory. So I just let them show up in the same directory where the file list was and then moved all the .gtfs

#default penalties for multimappers. 

# -u :Turn off multi-mapping correction. In the default case this correction is enabled, and each read that is mapped in n places only contributes 1/n to the transcript coverage instead of 1


