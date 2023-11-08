#Alriiiiiight. Let's see what we got here. 
#after filtering and trimming off the edges of SNP calls that aren't in the QTL region, we
#have 41000+ variants. Time to do some additional filtering with logical statements.

#We are interested in the following variants:

#Locations where...
# 02-08 = 02-52 = 03-08 && 04-12 = 02-19. Moreover, the first 3's genotypes should be 'contained' within the last two; i.e, the last two will have one bonus allele. 
# Balaton's genotype call should be contained within these individuals.
# and Surefire will have at least 1 matching 'character' or 'letter' with each of the lates. 
# don't love that because of such low coverage, 04-12 and 02-19 often have the same genotype but are genotyped incorrectly. i.e., 
# where balaton is G and surefire is C/G, their genotypes should be G/G/C; however, they usually get fucked up because G and C might
# have equal probabilities for the 3rd allele. So I'll get something shitty like C/C/G and G/G/G >.< We'll think about that more later. It might be a filter
# we apply much later than the others. 

#all of the early bloomers A subgenome alleles come from haplotype i EXCEPT FOR 03-01's in the first part of the QTL -- up until base pair 12608781. 
#then 03-01 has g and i. 

#i.e., we are interested in areas where 03-27 = 03-28 = 03-46 = 04-34 = 02-65

# | means the SNPs in the area are phased. We don't necessarily care about that right now.. supposedly they're all phased anyway.
# so I wonder if it would be better to just find and replace | with /
# I hesitate to drop places where there are missing values.. theoretically, if they have the same genotype, wouldn't a difficult area be difficult for all the individuals with that haplotype?
# I know shit happens but... we have 41000 variants. I should probably start out really being conservative. 

#now to code all this shit. 


# Side note: 04-12 and 02-19 differ in their sweet cherry genome composition, but have the same A haplotypes (k, k, and i); i is from surefire.

#Jan 30 2023: script updated after manual annotation of genes in the chromosome 4 QTL region. 
#And, updated again in March 2023 (round 2 update) since some of the DEGs that also had genetic variation looked a little funky. 

library(ggplot2)
library(tidyverse)
library(readr)

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/scripts_for_pop4_sequencing/final_scripts/subA_subA'_subB/variants/chr4_QTL_round2_update/")

filtered_variants = read.table("ALL_filtered_variants_4A_EFF_QTL_update2_table", header = TRUE)
summary(filtered_variants)

#first, let's drop some useless columns, because all of these variants passed our quality filters. If we want to know the quality of a certain variant later we can always unite the dataset by position.

filtered_variants2 = filtered_variants %>%
  dplyr::select(POS, TYPE, EFF, REF, ALT, BAL.GT, SF.GT, X27.02.08.GT, X27.02.52.GT, X27.03.08.GT, X27.02.19.GT, 
         X27.04.12.GT, X27.02.65.GT, X27.03.46.GT, X27.03.27.GT, X27.03.28.GT, X27.04.34.GT, X27.03.01.GT)

colnames(filtered_variants2) = c("POS", "TYPE", "PRED_EFFECT", "REF", "ALT", "balaton", "surefire", "02_08", "02_52", "03_08", "02_19", "04_12",
                                 "02_65", "03_46", "03_27", "03_28", "04_34", "03_01")

summary(filtered_variants2)
head(filtered_variants2)

#| is gonna mess with comparisons, so I'll cahnge it to / for now:
filtered_variants3 = data.frame(lapply(filtered_variants2, function(x) {
  gsub("\\|", "/", x)
}))

head(filtered_variants3) #sweet, it worked. But why is it adding Xs to my column headers...? To make sure we know it's a string? I'll leave it for now. 
summary(filtered_variants3)

#let's start applying logical statements.

#I feel ok with this statement; it should still put us closer to the causative variant(s) cuz it isn't as frequently fucked up like the 04-12 and 02-19 comparison.  

#these 3 late-bloomers all have 2X subgenome A, which are both from k. So, their calls should be identical.
filtered_variants4 = filtered_variants3[filtered_variants3$X02_08==filtered_variants3$X02_52 & filtered_variants3$X02_52==filtered_variants3$X03_08, ]

#all of the late-bloomers that have only i for their subgenome A and should have identical genotype calls / alleles. 
filtered_variants5 = filtered_variants4[filtered_variants4$X03_46==filtered_variants4$X03_28 &
                                          filtered_variants4$X03_28==filtered_variants4$X03_27 &
                                          filtered_variants4$X03_27==filtered_variants4$X02_65 &
                                          filtered_variants4$X02_65==filtered_variants4$X04_34, ]
head(filtered_variants5)

#a couple of late-bloomers are 3X for subgenome A; they have haplotypes k, k, and i. Therefore, they're genotype calls should be the same.
filtered_variants6B = filtered_variants5[filtered_variants5$X04_12==filtered_variants5$X02_19, ]

head(filtered_variants6B)
#while we're at it, we also don't want places where 02-19 and 04-12 are all X/X/X where X is the same allele -- this is because
#the genotype call / allele for i SHOULD BE DIFFERENT than k. 

#I am not sure how to code X where it is more than one character that's the same >.<
filtered_variants6C = filtered_variants6B %>%
  filter(X04_12 != "A/A/A") %>%
  filter(X04_12 != "T/T/T") %>%
  filter(X04_12 != "C/C/C") %>%
  filter(X04_12 != "G/G/G")

filtered_variants6C$POS = as.numeric(filtered_variants6C$POS)

head(filtered_variants6C)
summary(filtered_variants6C)

#furthermore... since we know no A subgenomes were inherited from Balaton for the early-bloomers, we can set those as unequal for further filtering. In other words, k will be different than the A subgenomes from i or g
#and we can compare them directly because they are the same ploidy (save for 03-01)
filtered_variants6D = filtered_variants6C[filtered_variants6C$balaton!=filtered_variants6C$X03_46, ]
ggplot(filtered_variants6D, aes(x=POS)) +
  geom_histogram()

#perhaps the thing to do now is to look for what kinds of variants these are with snpeff. Let's do that. 
#now, what kinds of variants are we working with? How many are Non-synonymous changes?

non_syn = dplyr::filter(filtered_variants6D, grepl("NON_SYNONYMOUS", PRED_EFFECT))

ggplot(non_syn, aes(x=POS)) +
  geom_histogram()
#still pretty scattered about. 

#any high effect variants?
high_effect = dplyr::filter(filtered_variants6D, grepl("HIGH", PRED_EFFECT))
ggplot(high_effect, aes(x=POS, y=0)) +
  geom_point(size=0.5)

ggplot(high_effect, aes(x=POS)) +
  geom_histogram()

#there are 13 of them. I think ENOD-like fell off the list when I reannotated it. 

modifying_effect = dplyr::filter(filtered_variants6D, grepl("MODIFIER", PRED_EFFECT))
#how many modifiers are classified as Downstream? Upstream? Intergenic? UTR? Intron?
nrow(modifying_effect[grepl("DOWNSTREAM", modifying_effect$PRED_EFFECT), ]) #1376 downstream variants
nrow(modifying_effect[grepl("UPSTREAM", modifying_effect$PRED_EFFECT), ]) #3801 upstream variants #perchance these are in the promoter region?
#nrow(modifying_effect[grepl("INTERGENIC", modifying_effect$PRED_EFFECT), ]) #902 intergenic varaints <- idk how we're going to functionally tie these to anything so imma ignore em for now. 
nrow(modifying_effect[grepl("UTR", modifying_effect$PRED_EFFECT), ]) #275 UTR variants
nrow(modifying_effect[grepl("INTRON", modifying_effect$PRED_EFFECT), ]) #259 intronic variants



write_excel_csv(non_syn, file="nonsynonymous_variants_of_k_haplotype_QTL_update2.csv")
write_excel_csv(high_effect, file="high_effect_variants_of_k_haplotype_QTL_update2.csv")
write_excel_csv(modifying_effect, file="modifier_variants_of_k_haplotype_QTL_update2.csv")






