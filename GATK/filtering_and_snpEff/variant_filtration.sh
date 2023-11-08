#!/bin/sh --login
#SBATCH -J filter
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64g
#SBATCH --time=1:00:00
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/combined_gvcf/variant_filter/variant_filter%j

cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/subset_4A/combined_gvcf/variant_filter

ml purge
ml GCCcore/10.2.0
ml GATK/4.2.0.0-Java-11
ml picard/2.25.0-Java-11

$EBROOTGATK/gatk --java-options "-Xmx60g" VariantFiltration \
-R /mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/Pcer_montAv5_p99_400k_scaffolded_assembly.fasta \
-V ../ALL_raw_variants_4A.vcf.gz \
-O ./ALL_filtered_variants_4A.vcf.gz \
--filter-name "QD_2" \
-filter "QD < 2.0" \
--filter-name "SOR_3" \
-filter "SOR > 3.0" \
--filter-name "MQ_50" \
-filter "MQ < 50.0" \
--filter-name "MQRankSum_neg8" \
-filter "DP < 2" \
--filter-name "DP_2" \
-filter "MQRankSum < -8.0" \
--filter-name "ReadPosRankSum_neg4" \
-filter "ReadPosRankSum < -4.0"


#don't worry about the error where the software doesn't recognize MQRankSum and ReadPosRankSum ---> https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/2017-01-18-2016-08-11/8323-Why-do-MQRankSum-and-ReadPosRankSum-not-appear-in-some-vcf-file-entries
#'Rank Sum Test annotations cannot be calculated for sites that are not heterozygous. They need to have a mix of ref and alt reads to be calculated'



#optional arguments of note:

#--missing-values-evaluate-as-failing
#When evaluating the JEXL expressions, missing values should be considered failing the expression... default is false. 

#--filter-expression / -filter
#One or more expression used with INFO fields to filter
#VariantFiltration accepts any number of JEXL expressions (so you can have two named filters by using --filterName One --filterExpression "X < 1" --filterName Two --filterExpression "X > 2").
#List[String]  []

#for notes on specs to hard-filter on:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants

#for general documentation of this tool:
# https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration

#quality by depth = QD
#This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. For filtering purposes it is better to use QD than either QUAL or DP directly.
#Hard filter recommendation: < 2

#RMSMappingQuality = MQ
#This is the root mean square mapping quality over all the reads at the site. Instead of the average mapping quality of the site, this annotation gives the square root of the average of the squares of the mapping qualities at the site. A low standard deviation means the values are all close to the mean, whereas a high standard deviation means the values are all far from the mean.When the mapping qualities are good at a site, the MQ will be around 60.
#Hard filter recommendation: < 40

#StrandOddsRatio = SOR 
#This is another way to estimate strand bias using a test similar to the symmetric odds ratio test. SOR was created because FS tends to penalize variants that occur at the ends of exons. Reads at the ends of exons tend to only be covered by reads in one direction and FS gives those variants a bad score. SOR will take into account the ratios of reads that cover both alleles.
#Hard filter recommendation: > 3
#Although there is a non-negligible population of variants with an SOR value less than 3 that failed VQSR, our hard filtering recommendation of failing variants with an SOR value greater than 3 will at least remove the long tail of variants that show fairly clear bias according to the SOR test. might wanna use this instead of FS to not overly penalize exons.

#FisherStrand = FS
#This is the Phred-scaled probability that there is strand bias at the site. Strand Bias tells us whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele. When there little to no strand bias at the site, the FS value will be close to 0. Might unfairly penalize reads mapping to the ends of exons

#MappingQualityRankSumTest = MQRankSum
#This is the u-based z-approximation from the Rank Sum Test for mapping qualities. It compares the mapping qualities of the reads supporting the reference allele and the alternate allele. A positive value means the mapping qualities of the reads supporting the alternate allele are higher than those supporting the reference allele; a negative value indicates the mapping qualities of the reference allele are higher than those supporting the alternate allele. A value close to zero is best and indicates little difference between the mapping qualities.
#Hard filter recommendation: [range]... not super clear. < -10.5 was mentioned. But what about the upper bound??

#ReadPosRankSumTest = ReadPosRankSum
#This is the u-based z-approximation from the Rank Sum Test for site position within reads. It compares whether the positions of the reference and alternate alleles are different within the reads. Seeing an allele only near the ends of reads is indicative of error, because that is where sequencers tend to make the most errors. A negative value indicates that the alternate allele is found at the ends of reads more often than the reference allele; a positive value indicates that the reference allele is found at the ends of reads more often than the alternate allele. A value close to zero is best because it indicates there is little difference between the positions of the reference and alternate alleles in the reads


#NOTE: Our hard filtering recommendations are meant to be very lenient. We prefer to keep all potentially decent variants rather than get rid of a few bad variants.
