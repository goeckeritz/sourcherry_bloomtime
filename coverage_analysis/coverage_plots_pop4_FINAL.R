#before we do any variant calling, we need to visualize what the hell is going on in terms of coverage for these pop4 individuals. 

#rm(list = ls(all.names = TRUE)) <- the oh shit button. Reset your environment and delete all objects.
setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/scripts_for_pop4_sequencing/final_scripts/subA_subA'_subB/pop4_coverage")

library(tidyverse)
library(ggplot2)
library(readr)

file_list = list.files(path = "./",
                       recursive = TRUE,
                       pattern = "\\.tsv$",
                       full.names = TRUE)
file_list


pop4_cov = readr::read_tsv(file_list, id = "tree")
pop4_cov$tree = gsub(".//", "", pop4_cov$tree)
pop4_cov$tree = gsub("_500000_headers.tsv", "", pop4_cov$tree)
pop4_cov$chr = gsub("Pcer_chr", "", pop4_cov$chr)
pop4_cov$chr = gsub("A__", "A'", pop4_cov$chr)
pop4_cov$chr = sub("([0-9])([A,B,A'])", "\\1_\\2", pop4_cov$chr, fixed=FALSE)


pop4_test = pop4_cov %>% separate(chr, c("chr_num", "homoeolog"), sep = "_")
pop4_test$homoeolog = as.factor(pop4_test$homoeolog)
pop4_test$chr_num = as.factor(pop4_test$chr_num)
pop4_test$tree = as.factor(pop4_test$tree)
levels(pop4_test$homoeolog)
levels(pop4_test$chr_num)
levels(pop4_test$tree)
str(pop4_test)
pop4_test=as.data.frame(pop4_test)

#cool! #now we can separate files by chr_num...

chr1 = dplyr::filter(pop4_test, chr_num == "1") %>% droplevels
chr2 = dplyr::filter(pop4_test, chr_num == "2") %>% droplevels
chr3 = dplyr::filter(pop4_test, chr_num == "3") %>% droplevels
chr4 = dplyr::filter(pop4_test, chr_num == "4") %>% droplevels
levels(chr4$chr_num)
chr5 = dplyr::filter(pop4_test, chr_num == "5") %>% droplevels
chr6 = dplyr::filter(pop4_test, chr_num == "6") %>% droplevels
chr7 = dplyr::filter(pop4_test, chr_num == "7") %>% droplevels
chr8 = dplyr::filter(pop4_test, chr_num == "8") %>% droplevels

no_parents = pop4_test %>%
  dplyr::filter(tree!="BAL") %>%
  dplyr::filter(tree!="SF") %>%
  droplevels()
levels(no_parents$tree)

parents = pop4_test[!grepl("27-0[2-4]", pop4_test$tree),] %>%
  droplevels()
levels(parents$tree) #fucking finally this is working. I cannot figure out why dplyr::filter was not working for this variable. 
str(no_parents)
str(parents)


chr1_no_parents = dplyr::filter(no_parents, chr_num == "1")%>% droplevels()
chr2_no_parents = dplyr::filter(no_parents, chr_num == "2")%>% droplevels()
chr3_no_parents = dplyr::filter(no_parents, chr_num == "3")%>% droplevels()
chr4_no_parents = dplyr::filter(no_parents, chr_num == "4") %>% droplevels()
chr5_no_parents = dplyr::filter(no_parents, chr_num == "5")%>% droplevels()
chr6_no_parents = dplyr::filter(no_parents, chr_num == "6")%>% droplevels()
chr7_no_parents = dplyr::filter(no_parents, chr_num == "7")%>% droplevels()
chr8_no_parents = dplyr::filter(no_parents, chr_num == "8")%>% droplevels()

chr1_parents = dplyr::filter(parents, chr_num == "1")%>% droplevels()
chr2_parents = dplyr::filter(parents, chr_num == "2")%>% droplevels()
chr3_parents = dplyr::filter(parents, chr_num == "3")%>% droplevels()
chr4_parents = dplyr::filter(parents, chr_num == "4")%>% droplevels()
chr5_parents = dplyr::filter(parents, chr_num == "5")%>% droplevels()
chr6_parents = dplyr::filter(parents, chr_num == "6")%>% droplevels()
chr7_parents = dplyr::filter(parents, chr_num == "7")%>% droplevels()
chr8_parents = dplyr::filter(parents, chr_num == "8")%>% droplevels()


#now we can look at coverage, one chromosome at a time.

layout.matrix=matrix(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), nrow = 3, ncol = 4)
layout.matrix

#where are those DAM genes.
ggplot(chr1_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog)) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,25) +
  geom_vline(xintercept = 43860000, linetype="dotted", 
             color = "blue4", size=0.5) + 
  geom_vline(xintercept = 43920000, linetype="dotted", 
             color = "blue4", size=0.5) +
  geom_vline(xintercept = 48700000, linetype="dotted", 
             color = "red3", size=0.5) + 
  geom_vline(xintercept = 48800000, linetype="dotted", 
             color = "red3", size=0.5) +
  geom_vline(xintercept = 46880000, linetype="dotted", 
             color = "green4", size=0.5) + 
  geom_vline(xintercept = 46990000, linetype="dotted", 
             color = "green4", size=0.5) 

ggplot(chr1_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,20) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 1") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))


ggplot(chr2_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,20) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 2") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))


ggplot(chr3_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,20) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 3") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

levels(chr4_no_parents$tree) #What the hell? why are Balaton and Surefire still there? If you look at the dataframe manually and SF and BAL aren't there anymore, like it's supposed to be. 
chr4_no_parents$tree = revalue(chr4_no_parents$tree, c("27-03-25"="Early 1","27-03-46"="Early 2","27-04-34"="Early 3","27-03-28"="Early 4","27-02-65"="Early 5","27-03-27" = "Early 6","27-03-01" = "Early 7","27-03-08"="Late 1","27-02-19"="Late 2","27-04-12"="Late 3","27-02-08" = "Late 4","27-02-52"="Late 6"))
chr4_no_parents$tree = factor(chr4_no_parents$tree, levels = c("Early 1","Early 2","Early 3","Early 4","Early 5","Early 6","Early 7","Late 1","Late 2","Late 3","Late 4","Late 6"))

ggplot(chr4_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), linewidth=0.75) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,20) +
  ylab("Coverage (mapped reads)") +
  xlab("Chromosome Position (Mb)") +
  scale_x_continuous(breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000), labels=c(0,5,10,15,20,25,30))+
  ggtitle("Chromosome 4") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick"), name="Subgenome") +
  theme_minimal() +
  theme(strip.text = element_text(size=12), axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=13), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=14, face="bold"), plot.title=element_text(size=20, face = "bold", hjust=0.5), legend.text=element_text(size=15, face="bold"), legend.title=element_text(size=16, face="bold", margin=margin(b=10))) +
  geom_vline(xintercept=8824991, color="black", linewidth=0.5, linetype="dotted")+
  geom_vline(xintercept=15886209, color="black", linewidth=0.5, linetype="dotted")

#each allele is supported by ~3-5 reads -- although the depth looks to be a little less in some individuals, so the dosage doesn't perfectly line up in a few cases. 
#e.g., Late 2, Early 2, Early 4, Early 6
ggplot(chr4_no_parents, aes(pos, cov/4)) +
  geom_line(aes(color=homoeolog), linewidth=0.75) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,4) +
  ylab("Relative Coverage") +
  xlab("Chromosome Position (Mb)") +
  scale_x_continuous(breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000), labels=c(0,5,10,15,20,25,30))+
  ggtitle("Chromosome 4") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick"), name="Subgenome") +
  theme_minimal() +
  theme(strip.text = element_text(size=12), axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=13), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=14, face="bold"), plot.title=element_text(size=20, face = "bold", hjust=0.5), legend.text=element_text(size=15, face="bold"), legend.title=element_text(size=16, face="bold", margin=margin(b=10))) +
  geom_vline(xintercept=8824991, color="black", linewidth=0.5, linetype="dotted")+
  geom_vline(xintercept=15886209, color="black", linewidth=0.5, linetype="dotted")

  

ggplot(chr5_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,20) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 5") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

ggplot(chr6_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,20) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 6") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

ggplot(chr7_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,20) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 7") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

ggplot(chr8_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,20) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 8") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))


ggplot(chr1_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog)) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,75) +
  geom_vline(xintercept = 43860000, linetype="dotted", 
             color = "blue4", size=0.5) + 
  geom_vline(xintercept = 43920000, linetype="dotted", 
             color = "blue4", size=0.5) +
  geom_vline(xintercept = 48700000, linetype="dotted", 
             color = "red3", size=0.5) + 
  geom_vline(xintercept = 48800000, linetype="dotted", 
             color = "red3", size=0.5) +
  geom_vline(xintercept = 46880000, linetype="dotted", 
             color = "green4", size=0.5) + 
  geom_vline(xintercept = 46990000, linetype="dotted", 
             color = "green4", size=0.5) 
  
  

ggplot(chr4_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog)) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,25) +
  geom_vline(xintercept = 8557351, linetype="dotted", 
           color = "blue4", size=1) + 
  geom_vline(xintercept = 15130451, linetype="dotted", 
             color = "blue4", size=1) +
  geom_vline(xintercept = 9276159, linetype="solid", 
           color = "green4", size=0.5) + 
  geom_vline(xintercept = 15886209, linetype="solid", 
             color = "green4", size=0.5) +
  geom_vline(xintercept = 8824991, linetype="dotted", 
             color = "red3", size=0.5) + 
  geom_vline(xintercept = 15130451, linetype="dotted", 
             color = "red3", size=0.5) +
  geom_vline(xintercept=13000000, linetype="solid",
             color = "transparent", size=0.2) +
  geom_vline(xintercept=12000000, linetype="solid",
           color = "transparent", size=0.2) +
  geom_vline(xintercept=12533308, linetype="solid",
             color="black", size=0.2)


#ggplot(chr4, aes(pos, cov)) +
 #geom_line(aes(color=homoeolog)) +
  #facet_wrap(vars(tree), nrow=4, ncol=4) +
  #ylim(0,25)

ggplot(chr4_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog)) +
  facet_wrap(vars(tree), nrow=4, ncol=4) +
  ylim(0,75) +
  geom_vline(xintercept = 8557351, linetype="dotted", 
             color = "blue4", size=1) + 
  geom_vline(xintercept = 15130451, linetype="dotted", 
             color = "blue4", size=1) +
  geom_vline(xintercept = 9276159, linetype="solid", 
             color = "green4", size=0.5) + 
  geom_vline(xintercept = 15886209, linetype="solid", 
             color = "green4", size=0.5) +
  geom_vline(xintercept = 8824991, linetype="dotted", 
             color = "red3", size=0.5) + 
  geom_vline(xintercept = 15130451, linetype="dotted", 
             color = "red3", size=0.5) +
  geom_vline(xintercept=12533308, linetype="solid",
             color="black", size=0.2)

#colocalizing marker for QTL in pop4:
#RosBREED_snp_tart_cherry_f_Pp4_10832168
#flanking markers found on the genetic map: https://www.rosaceae.org/mapviewer/370/LG4
#RosBREED_snp_tart_cherry_f_Pp4_07147868
#RosBREED_snp_tart_cherry_a_Pp4_12532690

#chr_4B
#colocalizing: 13555134 (best match, but only 92.6% identity)
#8557351, 8562280 (maps to two places...)
#15130451

#chr_4A:
#colocalizing: 13840124 (best match, but only 90.2% identity) <- that is REALLY close to that ENOD-like gene. it's less than 200 kb away
#8824991
#15886209
  
#chr_4A':
#colocalizing:	14001951 (best match, but only 95.6% identity)
#9276159
#15705568

#colocalizing markers had a alignment length of 184-195

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

find_split = read.table("../subset_27-03-01_chr4A.txt", header=TRUE)

ggplot(find_split, aes(pos, cov)) +
  geom_line(color="hotpink") +
  ylim(0,25) +
  geom_vline(xintercept=12500000, linetype="solid",
              color="black", size=0.2) +
  geom_vline(xintercept=12600000, linetype="solid",
             color="black", size=0.2) +
  geom_vline(xintercept=12533308, linetype="dashed",
             color="darkgreen", size=0.2)

#split is very close to 12.5 Mb, just a liiiiiitle more
#split is just about at 12533308 :)


ggplot(chr4_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1) +
  facet_wrap(vars(tree), nrow=4, ncol=4) +
  ylim(0,75) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 4") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  geom_vline(xintercept = 8824991, linetype="dotted", 
             color = "black", size=1) + 
  geom_vline(xintercept = 15130451, linetype="dotted", 
             color = "black", size=1) +
  theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))



ggplot(chr4_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1.5) +
  facet_wrap(vars(tree), nrow=4, ncol=4) +
  ylim(0,75) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 4") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  geom_vline(xintercept = 8824991, linetype="dotted", 
             color = "black", size=1) + 
  geom_vline(xintercept = 15130451, linetype="dotted", 
             color = "black", size=1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))


ggplot(chr4_no_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1) +
  facet_wrap(vars(tree), nrow=3, ncol=4) +
  ylim(0,20) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 4") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  geom_vline(xintercept = 8824991, linetype="dotted", 
             color = "black", size=1) + 
  geom_vline(xintercept = 15130451, linetype="dotted", 
             color = "black", size=1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))


ggplot(chr6_parents, aes(pos, cov)) +
  geom_line(aes(color=homoeolog), size=1.5) +
  facet_wrap(vars(tree), nrow=4, ncol=4) +
  ylim(0,75) +
  ylab("Coverage") +
  xlab("Chromosome Position") +
  ggtitle("Chromosome 6") +
  scale_color_manual(values=c("skyblue", "palegreen", "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))
