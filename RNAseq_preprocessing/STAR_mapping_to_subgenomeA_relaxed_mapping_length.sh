#!/bin/sh --login
#SBATCH -J STARXS_subgenomeA_relaxed
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=12g
#SBATCH --time=48:00:00
#SBATCH -o /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/STAR/alignment_subgenomeA_relaxed_%j
#SBATCH -a 1-48
#SBATCH --export=INFILE1=/mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/trimmomatic/paired_reads/RNAseq_paired.txt

module purge
module use /mnt/home/johnj/software/modulefiles
module load GCC/7.3.0-2.30  OpenMPI/3.1.1-CUDA
module load STAR/2.7.3a

subgenomeA=/mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/chr4_QTL_update/subgenomeA.fasta
subgenomeA_index=/mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/STAR_index_subA/
#gtf_scaffolded= <-- used when building the STAR index to make splice junction alignments more accurate. Highly recommended in the documentation. 

FILE1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE1}`
echo ${FILE1}
Read1_compressed=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE1}`
echo ${Read1_compressed}
Read2_compressed=$(echo ${Read1_compressed}| sed "s/__R1/__R2/")
echo ${Read2_compressed}
Name=$(echo ${FILE1}| sed 's/\/mnt\/scratch\/goeckeri\/RNAseq_bloom\/20230210_mRNASeq_PE150\/trimmomatic\/paired_reads\///' | sed "s/__R1_.*//")
echo ${Name}
Read1=$(echo ${Read1_compressed}| sed "s/\.gz//")
echo ${Read1}
Read2=$(echo ${Read2_compressed}| sed "s/\.gz//")
echo ${Read2}

cd /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/trimmomatic/paired_reads/
gunzip -c ${Read1_compressed} > ${Read1}
gunzip -c ${Read2_compressed} > ${Read2}


cd /mnt/scratch/goeckeri/RNAseq_bloom/20230210_mRNASeq_PE150/STAR/subgenomeA/subgnomeA_relaxed/
STAR --runMode alignReads --runThreadN 28 --genomeDir $subgenomeA_index --readFilesIn ${Read1} ${Read2} --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 10 --outSAMattributes NH HI AS nM XS --outFilterMismatchNoverLmax 0.3 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --alignIntronMax 7000 --alignSoftClipAtReferenceEnds Yes --outFileNamePrefix ${Name}_subA_relaxed --outTmpDir ${Name}_tmp --outReadsUnmapped Fastx --outStd BAM_SortedByCoordinate > ${Name}_subA_relaxed_sorted.bam

# | samtools view -@ 28 -Su | samtools sort -@ 28 -O bam -o ${Name}_aln.sorted.bam
#samtools index ${Name}_aln.sorted.bam

#I thought about adding the --outSAMprimaryFlag to this run, and changing it from the default to AllBestScore, thinking it would be good to give alignments of a alternative allele a fighting chance to also be called 'primary.' However, there's only going to be one of them in the reference.. and I wouldn't want to confuse future software like stringtie and cufflinks and CLASS2 when making the transcript model. 
#So, it would thus be best to represent one allele as accurately as we can.
#--outSAMmultNmax == max number of multiple alignments for a read that will be output to the SAM/BAM files. 10 is plenty in my opinion. 
#the options for removing duplicates was a bit weird, so I probably should use picard tools or something for that. 
#For the sake of ploidy and limiting multi-mappers, I applied --outFilterMismatchNoverLmax. Documentation says: alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value. 
#--outFilterScoreMinOverLread = alignment will only be output if it has a higher score than this value, which is normalized to the read length (sum of mates' lengths for paired-end reads). I upped this to 0.75 to make the output more conservative -- default is 0.66. 
#I should investigate more what penalizes a mapping score.. mismatches, how much clipping, and multimapping, probably. see --score* options. 
#alignIntronMax = maximum size of expected introns. 
# --readStrand was a mistake in the manual. STAR doesn't even utilize it... https://github.com/alexdobin/STAR/issues/818

#when I aligned 02-19 to the I got like 30% unmapped reads due to the reason of 'too short' (this means the aligned length was too short); I'm not the only one whose ever had this problem of course. https://github.com/alexdobin/STAR/issues/169
#basically bostanict's fix was to do filtering of really shitty alignments with samtools or something. Maybe I'll end up filtering these out in the end anyway >.<
#Really, I should pay attention to average mapped length to see how much it goes down when these parameters are relaxed. 

#see also: https://github.com/alexdobin/STAR/issues/164

#I'm controlling mismatches by using outFilterMismatchNoverLmax -- I put it back at the default, 0.3; meaning whatever the length aligned, only 30% mismatches are allowed. 

