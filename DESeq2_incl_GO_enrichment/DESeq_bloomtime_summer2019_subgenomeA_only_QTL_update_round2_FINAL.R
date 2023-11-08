######Doing gene-level analyses for early and late bloomers#######
#syntelog lists created with GENESPACE - seemed to get a few thousand more as compared with MCScan
#(but double check this of course, the results still require a bit of interpretation)
#this script utilizes the second round of updates of the annotation I did for the QTL region. 

#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tibble)
library(pheatmap)
library(seriation)
library(dendextend)

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/DESeq2/QTL_round2_update/")
subA_transcript_countData <- as.matrix(read.csv("transcript_count_matrix.csv", row.names="transcript_id"))

colnames(subA_transcript_countData) <- sub("^X","",perl = T, colnames(subA_transcript_countData))
colnames(subA_transcript_countData) <- gsub(".","_",fixed = TRUE, colnames(subA_transcript_countData))
head(subA_transcript_countData)

nrow(subA_transcript_countData)
ncol(subA_transcript_countData)

subA_countsNonZero <- subA_transcript_countData[apply(subA_transcript_countData,1,function(x){!all(x == 0)}),]
nrow(subA_countsNonZero)

empty_rows = data.frame(matrix(ncol=ncol(subA_countsNonZero), nrow=3))
colnames(empty_rows) = colnames(subA_countsNonZero)
subA_w_factors = rbind(subA_countsNonZero, empty_rows)

for (i in 1:ncol(subA_w_factors)){
  if (grepl("27_02_19",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_02_19"
  } else if (grepl("27_03_08",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_03_08"
  } else if (grepl("27_03_25",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_03_25"
  } else if (grepl("27_03_46",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_03_46"
  } else if (grepl("27_04_12",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_04_12"
  } else if (grepl("27_04_34",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_04_34"
  }
  if (grepl("6_12_19",colnames(subA_w_factors)[i])){
    subA_w_factors[nrow(subA_w_factors)-1,i] <- "Jun12"
  } else if (grepl("6_29_19",colnames(subA_w_factors)[i])){
    subA_w_factors[nrow(subA_w_factors)-1,i] <- "Jun29"
  }  else if (grepl("7_18_19",colnames(subA_w_factors)[i])){
    subA_w_factors[nrow(subA_w_factors)-1,i] <- "Jul18"
  }
  if (grepl("27_02_19",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "late"
  } else if (grepl("27_03_08",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "late"
  } else if (grepl("27_03_25",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "early"
  } else if (grepl("27_03_46",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "early"
  } else if (grepl("27_04_12",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "late"
  } else if (grepl("27_04_34",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "early"
  }
}
  
row.names(subA_w_factors)[(nrow(subA_w_factors)-2):(nrow(subA_w_factors))] <- c("tree","date","bloom")
subA_coldata <- subA_w_factors[(nrow(subA_w_factors)-2):(nrow(subA_w_factors)),]
subA_coldata <- data.frame(t(subA_coldata))
all(rownames(subA_coldata) == colnames(subA_countsNonZero))
  
  #to prep individual tree to be nested in our linear model. 
for (i in 1:nrow(subA_coldata)){
    if (grepl("27_02_19", subA_coldata$tree[i])){
      subA_coldata[i, "nested_individual"] <- "1"
    } else if (grepl("27_03_08", subA_coldata$tree[i])){
      subA_coldata[i, "nested_individual"] <- "2"
    } else if (grepl("27_04_12", subA_coldata$tree[i])){
      subA_coldata[i, "nested_individual"] <- "3"
    } else if (grepl("27_03_25", subA_coldata$tree[i])){
      subA_coldata[i, "nested_individual"] <- "1"
    } else if (grepl("27_03_46", subA_coldata$tree[i])){
      subA_coldata[i, "nested_individual"] <- "2"
    } else if (grepl("27_04_34", subA_coldata$tree[i])){
      subA_coldata[i, "nested_individual"] <- "3"
    }
  }
  
subA_coldata$nested_individual = as.factor(subA_coldata$nested_individual)
  
subA_dds3 <- DESeqDataSetFromMatrix(subA_countsNonZero, 
                                      colData=subA_coldata, 
                                      design= ~ bloom + bloom:nested_individual + bloom:date + date)
#The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input.  
#https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html <-- DESeq2 does median of ratios. Assumes that gene lengths are constant across samples, and they're all against the same genome so.. that should work. 
#https://www.reneshbedre.com/blog/expression_units.html#deseq-or-deseq2-normalization-median-of-ratios-method
#https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html
vsdC <- rlog(subA_dds3)
norm.countsC <- assay(vsdC)
head(norm.countsC[,1:6])
counts.corC <- cor(norm.countsC,method = "pearson")
  
pheatmap(counts.corC) #lol dat outlier. I wonder why it doesn't show up when mapping to any other subgenomes. 
  
baseline_PCAc = plotPCA(vsdC,intgroup=c("bloom", "nested_individual", "date"), returnData=TRUE)
percentVarC <- round(100 * attr(baseline_PCAc, "percentVar"))
plotPCA(vsdC,intgroup=c("bloom", "nested_individual", "date"))
head(baseline_PCAc)
  
levels(baseline_PCAc$group)
baseline_PCAc$for_plotting = NA
for (i in 1:nrow(baseline_PCAc)){
  if (grepl("early:1", baseline_PCAc$group[i])){
    baseline_PCAc[i, "for_plotting"] <- "Early 1"
  } else if (grepl("early:2", baseline_PCAc$group[i])){
    baseline_PCAc[i, "for_plotting"] <- "Early 2"
  } else if (grepl("early:3", baseline_PCAc$group[i])){
    baseline_PCAc[i, "for_plotting"] <- "Early 3"
  } else if (grepl("late:1", baseline_PCAc$group[i])){
    baseline_PCAc[i, "for_plotting"] <- "Late 1"
  } else if (grepl("late:2", baseline_PCAc$group[i])){
    baseline_PCAc[i, "for_plotting"] <- "Late 2"
  } else if (grepl("late:3", baseline_PCAc$group[i])){
    baseline_PCAc[i, "for_plotting"] <- "Late 3"
  }
}
baseline_PCAc$for_plotting = as.factor(baseline_PCAc$for_plotting)
baseline_PCAc$date = as.factor(baseline_PCAc$date)
baseline_PCAc$date = factor(baseline_PCAc$date, levels = c('Jun12','Jun29','Jul18'))

ggplot(baseline_PCAc, aes(PC1, PC2, color=for_plotting, shape=for_plotting, size=date)) +
  geom_point(stroke=1) +
  ggtitle("PCA for reads mapped to subgenome A")+
  xlab(paste0("PC1: ",percentVarC[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarC[2],"% variance")) + 
  scale_color_manual(values=c("steelblue","steelblue2","steelblue1","pink2","violetred2","violetred3"), name="Individual")+
  scale_shape_manual(values=c(1,2,5,1,2,5), name="Individual", labels=c("Early 1","Early 2","Early 3","Late 1","Late 2", "Late 3")) +
  scale_size_manual(values=c(3, 5, 9), name="Date", labels=c("Jun12"="June 12", "Jun29"="June 29", "Jul18"="July 18")) +
  theme_minimal()+
  theme(axis.text.x = element_text(hjust = 1, size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15), 
        plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold", margin=margin(b=10))) 
  
#####let's drop the potential outliers. 
  
subA_countsNonZero_df <- as.data.frame(subA_countsNonZero)
subA_countsNonZero_no_outlier <- as.matrix(subA_countsNonZero_df[,!names(subA_countsNonZero_df) %in% c("27_04_34_6_12_19_R2","27_03_46_6_29_19_R3")])
head(subA_countsNonZero_no_outlier)
  
subA_coldata_no_outlier <- subA_coldata[!row.names(subA_coldata) %in% c("27_04_34_6_12_19_R2","27_03_46_6_29_19_R3"),]
all(rownames(subA_coldata_no_outlier) == colnames(subA_countsNonZero_no_outlier))
  
#rerun the model once the outliers are dropped. 
subA_dds3b <- DESeqDataSetFromMatrix(subA_countsNonZero_no_outlier, 
                                       colData=subA_coldata_no_outlier, 
                                       design= ~ bloom + bloom:nested_individual + bloom:date + date)
  
vsdX <- rlog(subA_dds3b)
norm.countsX <- assay(vsdX)
head(norm.countsX[,1:6])
counts.corX <- cor(norm.countsX,method = "pearson")
head(counts.corX)
pheatmap(counts.corX)
 
baseline_PCAX = plotPCA(vsdX,intgroup=c("bloom", "nested_individual", "date"), returnData=TRUE)
percentVarX <- round(100 * attr(baseline_PCAX, "percentVar"))
plotPCA(vsdX,intgroup=c("bloom", "nested_individual", "date"))
head(baseline_PCAX)
  
levels(baseline_PCAX$group)
baseline_PCAX$for_plotting = NA
for (i in 1:nrow(baseline_PCAX)){
  if (grepl("early:1", baseline_PCAX$group[i])){
    baseline_PCAX[i, "for_plotting"] <- "Early 1"
  } else if (grepl("early:2", baseline_PCAX$group[i])){
    baseline_PCAX[i, "for_plotting"] <- "Early 2"
  } else if (grepl("early:3", baseline_PCAX$group[i])){
    baseline_PCAX[i, "for_plotting"] <- "Early 3"
  } else if (grepl("late:1", baseline_PCAX$group[i])){
    baseline_PCAX[i, "for_plotting"] <- "Late 1"
  } else if (grepl("late:2", baseline_PCAX$group[i])){
    baseline_PCAX[i, "for_plotting"] <- "Late 2"
  } else if (grepl("late:3", baseline_PCAX$group[i])){
    baseline_PCAX[i, "for_plotting"] <- "Late 3"
  }
}
baseline_PCAX$for_plotting = as.factor(baseline_PCAX$for_plotting)
baseline_PCAX$date = as.factor(baseline_PCAX$date)
baseline_PCAX$date = factor(baseline_PCAX$date, levels = c('Jun12','Jun29','Jul18'))

ggplot(baseline_PCAX, aes(PC1, PC2, color=for_plotting, shape=for_plotting, size=date)) +
  geom_point(stroke=1) +
  ggtitle("PCA for reads mapped to subgenome A")+
  xlab(paste0("PC1: ",percentVarX[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarX[2],"% variance")) + 
  scale_color_manual(values=c("steelblue","steelblue2","steelblue1","pink2","violetred2","violetred3"), name="Individual")+
  scale_shape_manual(values=c(1,2,5,1,2,5), name="Individual", labels=c("Early 1","Early 2","Early 3","Late 1","Late 2", "Late 3")) +
  scale_size_manual(values=c(3, 5, 9), name="Date", labels=c("Jun12"="June 12", "Jun29"="June 29", "Jul18"="July 18")) +
  theme_minimal()+
  theme(axis.text.x = element_text(hjust = 1, size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15), 
        plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold", margin=margin(b=10))) 

#ehhhhh.. 27-03-46_6_29_19_R3 came up as an outlier in the WGCNA analysis. Looking at it here again, 
#it might be good to drop it here too. DONE
  
#Side adventure: a heatmap of genes expressed in the QTL region. No DE analyses done yet. 
#write.table(subA_countsNonZero_no_outlier,"subA_no_outlier_counts_non_zero.table")

#subA_subset = read.csv("QTL_subset_count_matrix.csv", row.names="transcript_id")
#need to drop the outlier's column
#subA_subset = subA_subset[, colnames(subA_subset)[colnames(subA_subset) != 'X27.04.34_6_12_19_R2']] 
#subA_subset = subA_subset[, colnames(subA_subset)[colnames(subA_subset) != 'X27.03.46_6_29_19_R3']] 
 

#head(subA_subset)
#summary(subA_subset)
#colnames(subA_subset) <- sub("^X","",perl = T, colnames(subA_subset))
#colnames(subA_subset) <- gsub(".","_",fixed = TRUE, colnames(subA_subset))
#pheatmap(subA_subset) #woof, maybe we normalize the counts first. 

#subA_subset_dds <- DESeqDataSetFromMatrix(subA_subset, 
                                          #colData=subA_coldata_no_outlier, 
                                          #design= ~ bloom + bloom:nested_individual + bloom:date + date)

#vsdA = rlog(subA_subset_dds)
#norm.countsA = assay(vsdA)
#pheatmap(norm.countsA) #oof, need smaller chunks to look at.
#head(norm.countsA)


#p1_norm.countsA = norm.countsA[1:100,]
#p2_norm.countsA = norm.countsA[101:201,]
#p3_norm.countsA = norm.countsA[202:302,]
#p4_norm.countsA = norm.countsA[303:403,]
#p5_norm.countsA = norm.countsA[404:504,]
#p6_norm.countsA = norm.countsA[505:601,]
#pheatmap(p1_norm.countsA) 
#pheatmap(p2_norm.countsA) #not seeing an obvious pattern for any of these genes between earlies and lates. 
#pheatmap(p3_norm.countsA) #maybe Pcer_024813?
#pheatmap(p4_norm.countsA) #not seeing a discernible pattern. 
#pheatmap(p5_norm.countsA) #looks mostly like a date separation. Maybe Pcer_024910 showing something? but not at all dates. 
#pheatmap(p6_norm.countsA) #Pcer_024709? very subtle 
#sanity check. We should get a grouping of earlies in lates in just a PCA of these genes. 
#baseline_PCAA = plotPCA(vsdA,intgroup=c("bloom", "nested_individual", "date"), returnData=TRUE)
#percentVarA <- round(100 * attr(baseline_PCAA, "percentVar"))
#ggplot(baseline_PCAA, aes(PC1, PC2, color=bloom, shape=nested_individual, size=date)) +
#  geom_point(stroke=1.5) +
#  ggtitle("QTL region, mapped to subgenomeA only") +
#  xlab(paste0("PC1: ",percentVarA[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVarA[2],"% variance")) + 
#  scale_color_manual(values=c("late"="pink", "early"="skyblue")) +
#  scale_shape_manual(values=c(1,2,5)) +
#  scale_size_manual(values=c(3, 9, 5)) +
#  theme_minimal()


# DESeq2: Calculate the test statistics -----------------------------------

subA_dds = estimateSizeFactors(subA_dds3b)
subA_dds = estimateDispersions(subA_dds)
subA_dds = nbinomWaldTest(subA_dds, maxit=100000)

plotDispEsts(subA_dds) 
resultsNames(subA_dds) #da hell are these factors even. They seem to be a random bunch of comparisons 
#we need to relevel our factors. Idk why they end up like this. 

subA_dds$bloom <- relevel(subA_dds$bloom,ref = "late") #redefines what the 'control' is, and defines it as the reference comparison 
subA_dds$date <- relevel(subA_dds$date,ref="Jun12") #so our baselevel will be the late-bloomers on Jun 12. 
subA_dds$group <- factor(paste0(subA_dds$bloom,subA_dds$date)) #paste 0 is just concatenating the strings. We're going to be able to 
#look at the bloom groups in the individual bloom groups now. Reshaping the design with grouped variables:
design(subA_dds) <- ~group
head(subA_dds)
#I want to know what changes between jun12 and jun29, and june29 and july18. I could do june12 v july18 but would that really be all that informative? 
#Not sure. And if I do the first two contrasts, the third isn't orthagonal. Really the PCA for subA indicates
#the bloom groups are different at each date. So given the same environment, maybe k is modulating the 
#bloom time all of the time.

#since we regrouped the variables, we need to estimate the variance again:
subA_dds = DESeq(subA_dds) ###DESeq object for GO enrichment
plotDispEsts(subA_dds)
str(subA_dds)
subA_dds$group

# DESeq2: Extract and Filter the Results ---------------------------------- Contrast time, bebe!

early_v_late_Jun12 <- results(subA_dds,contrast = c("group","earlyJun12","lateJun12"),alpha = 0.05,format = "DataFrame") 
early_v_late_Jun29 <- results(subA_dds,contrast = c("group","earlyJun29","lateJun29"),alpha = 0.05,format = "DataFrame") 
early_v_late_Jul18 <- results(subA_dds,contrast = c("group","earlyJul18","lateJul18"),alpha = 0.05,format = "DataFrame") 
late_Jun12_Jun29 <- results(subA_dds,contrast = c("group","lateJun29","lateJun12"),alpha = 0.05,format = "DataFrame") 
late_Jun29_Jul18 <- results(subA_dds,contrast = c("group","lateJul18","lateJun29"),alpha = 0.05,format = "DataFrame") 
early_Jun12_Jun29 <- results(subA_dds,contrast = c("group","earlyJun29","earlyJun12"),alpha = 0.05,format = "DataFrame") 
early_Jun29_Jul18 <- results(subA_dds,contrast = c("group","earlyJul18","earlyJun29"),alpha = 0.05,format = "DataFrame") 


dds.res <- list(
  "early_v_late_Jun12"=early_v_late_Jun12,
  "early_v_late_Jun29"=early_v_late_Jun29,
  "early_v_late_Jul18"=early_v_late_Jul18,
  "late_Jun12_Jun29"=late_Jun12_Jun29,
  "late_Jun29_Jul18"=late_Jun29_Jul18,
  "early_Jun12_Jun29"=early_Jun12_Jun29,
  "early_Jun29_Jul18"=early_Jun29_Jul18
)

#building a little filtering function weeeeeee; setting the log fold change more relaxed as the gene(s) responsible for k may be quite subtly different. 
dds.filter <- function(x){
  x <- x[!is.na(x$padj),]
  y <- x[((x$padj < 0.05) & (abs(x$log2FoldChange) > 0.5)),]
  return(y)
}

#now apply the filter to all the dataframes in our list. 
dds.res.flt <- lapply(dds.res,dds.filter)
head(dds.res.flt)
lapply(dds.res.flt, summary)
#all my contrasts are here... 

# DESeq2: Plot the Log Fold Changes Against the Mean Counts ---------------
par(mfrow=c(2,4))
drawlines <- function() abline(h=c(-0.5,0.5),col="pink",lwd=2)

plotMA(dds.res$early_v_late_Jun12, alpha=0.05,main="Early- vs Late-bloomers on June 12th"); drawlines()
plotMA(dds.res$early_v_late_Jun29, alpha=0.05,main="Early- vs Late-bloomers on June 29th"); drawlines()
plotMA(dds.res$early_v_late_Jul18, alpha=0.05,main="Early- vs Late-bloomers on July 18th"); drawlines()
plotMA(dds.res$late_Jun12_Jun29, alpha=0.05,main="Changes in the late-bloomers from June 12th to June 29th"); drawlines()
plotMA(dds.res$early_Jun12_Jun29, alpha=0.05,main="Changes in the early-bloomers from June 12th to June 29th"); drawlines()
plotMA(dds.res$late_Jun29_Jul18, alpha=0.05,main="Changes in the late-bloomers from June 29th to July 18th"); drawlines()
plotMA(dds.res$early_Jun29_Jul18, alpha=0.05,main="Changes in the early-bloomers from June 29th to July 18th"); drawlines()

par(mfrow=c(1,1))
#mmmkay..

funct_anno = read.table("functional_annotation_QTL_round2_REFORMATTED_B.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE)

funct_anno = funct_anno %>%
  dplyr::select(Gene_ID,"Putative_Function"="PANTHER...Pfam")
head(funct_anno)

#anyway, back to filtering. 
for (i in 1:length(dds.res.flt)){
  dds.res.flt[[i]]$Gene_ID <- row.names(dds.res.flt[[i]])
}


###can add a direction down below for GO enrichment later

head(dds.res.flt$early_v_late_Jun12)

dds.res.flt2 <- lapply(dds.res.flt,function(x,y){left_join(as.data.frame(x),y,by="Gene_ID")},y=funct_anno)
head(dds.res.flt2$early_v_late_Jun12)
str(dds.res.flt2)  

write_excel_csv(dds.res.flt2$early_v_late_Jun12, "early_v_late_Jun12_DEGs.csv")
write_excel_csv(dds.res.flt2$early_v_late_Jun29, "early_v_late_Jun29_DEGs.csv")
write_excel_csv(dds.res.flt2$early_v_late_Jul18, "early_v_late_Jul18_DEGs.csv")
write_excel_csv(dds.res.flt2$late_Jun12_Jun29, "late_Jun12_Jun29_DEGs.csv")
write_excel_csv(dds.res.flt2$early_Jun12_Jun29, "early_Jun12_Jun29_DEGs.csv")
write_excel_csv(dds.res.flt2$late_Jun29_Jul18, "late_Jun29_Jul18_DEGs.csv")
write_excel_csv(dds.res.flt2$early_Jun29_Jul18, "early_Jun29_Jul18_DEGs.csv")


#let's just look at the QTL region genes now. 
QTL_genes <- read.table("QTL_gene_list.txt", header=TRUE)
head(QTL_genes)
xxx <- right_join(funct_anno, QTL_genes, by="Gene_ID") #keeps all QTL genes in this list, even if they have NA for their functional domains. 
nrow(xxx)
head(xxx)

dds.res.flt_QTL4 <- lapply(dds.res.flt2,function(x,y){inner_join(as.data.frame(x),y,by="Gene_ID")},y=xxx) #overlaps our full genome DEGs with genes in the QTL region, attaching functional information where available. 
head(dds.res.flt_QTL4$early_v_late_Jun12) 
lapply(dds.res.flt_QTL4, summary)

write_excel_csv(dds.res.flt_QTL4$early_v_late_Jun12, "early_v_late_Jun12_DEGs_chr4_QTL.csv")
write_excel_csv(dds.res.flt_QTL4$early_v_late_Jun29, "early_v_late_Jun29_DEGs_chr4_QTL.csv")
write_excel_csv(dds.res.flt_QTL4$early_v_late_Jul18, "early_v_late_Jul18_DEGs_chr4_QTL.csv")
write_excel_csv(dds.res.flt_QTL4$late_Jun12_Jun29, "late_Jun12_Jun29_DEGs_chr4_QTL.csv")
write_excel_csv(dds.res.flt_QTL4$early_Jun12_Jun29, "early_Jun12_Jun29_DEGs_chr4_QTL.csv")
write_excel_csv(dds.res.flt_QTL4$late_Jun29_Jul18, "late_Jun29_Jul18_DEGs_chr4_QTL.csv")
write_excel_csv(dds.res.flt_QTL4$early_Jun29_Jul18, "early_Jun29_Jul18_DEGs_chr4_QTL.csv")

list2env(dds.res.flt_QTL4,envir=.GlobalEnv) #NOTE - this splits the QTL DEGs list into its corresponding dataframes.
#if you want to do this for the full DEG list of the whole genome, change the object to dds.res.flt2, and note that the dataframes will be
#overwritten and will now be the DEG lists for the whole genome. 


###for heatmap visualization:
dev.off()

z <- full_join(early_v_late_Jun12,early_v_late_Jun29, by="Gene_ID", suffix=c("_early_v_late_Jun12", "_early_v_late_Jun29")) %>%
  dplyr::select("Gene_ID", "log2FoldChange_early_v_late_Jun12", "log2FoldChange_early_v_late_Jun29")
head(z)
zz <- early_v_late_Jul18 %>%
  dplyr::select("Gene_ID", "log2FoldChange_early_v_late_Jul18"="log2FoldChange")



DEGS_bloom_per_date <- full_join(z,zz, by="Gene_ID")
head(DEGS_bloom_per_date)
colnames(DEGS_bloom_per_date) = c("Gene_ID","June 12","June 29","July 18")
write_csv(DEGS_bloom_per_date, "DEGs_any_bloom_date.csv")

head(DEGS_bloom_per_date)
June12_ggheatmap = DEGS_bloom_per_date %>%
  dplyr::select("Gene_ID","log2FoldChange"="June 12") %>%
  drop_na()
  
June12_ggheatmap$Gene2 <- gsub("-R[A,B]","",June12_ggheatmap$Gene_ID)
June12_ggheatmap$Gene2 <- reorder(June12_ggheatmap$Gene2, June12_ggheatmap$log2FoldChange)
head(June12_ggheatmap)

ggplot(June12_ggheatmap, aes(x=0, Gene2, fill=log2FoldChange)) +
  geom_tile(aes(width=10), color="black") +
  ylab("Gene ID")+
  xlab("")+
  scale_fill_gradient2(low="cornflowerblue", mid="white", high="firebrick", name="", limits=c(-10,10), labels=c())+
    #guide = guide_legend(label.position="bottom")) +
  coord_flip()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=90, face="bold"), axis.title.x=element_text(size=15, angle=180),
        axis.text.y=element_text(size=0), legend.position = "top", legend.key.width=unit(4,'cm'), legend.text=element_text(size=12, angle=90))

head(DEGS_bloom_per_date)
June29_ggheatmap = DEGS_bloom_per_date %>%
  dplyr::select("Gene_ID","log2FoldChange"="June 29") %>%
  drop_na()

June29_ggheatmap$Gene2 <- gsub("-R[A,B]","",June29_ggheatmap$Gene_ID)
June29_ggheatmap$Gene2 <- reorder(June29_ggheatmap$Gene2, June29_ggheatmap$log2FoldChange)


ggplot(June29_ggheatmap, aes(x=0, Gene2, fill=log2FoldChange)) +
  geom_tile(aes(width=10), color="black") +
  ylab("Gene ID")+
  xlab("")+
  scale_fill_gradient2(low="cornflowerblue", mid="white", high="firebrick", name="", limits=c(-10,10), labels=c())+
  #guide = guide_legend(label.position="bottom")) +
  coord_flip()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=90, face="bold", size=12), axis.title.x=element_text(size=15, angle=180),
        axis.text.y=element_text(size=0), legend.position = "top", legend.key.width=unit(4,'cm'), legend.text=element_text(size=12, angle=90))

head(DEGS_bloom_per_date)
July18_ggheatmap = DEGS_bloom_per_date %>%
  dplyr::select("Gene_ID","log2FoldChange"="July 18") %>%
  drop_na()

July18_ggheatmap$Gene2 <- gsub("-R[A,B]","",July18_ggheatmap$Gene_ID)
July18_ggheatmap$Gene2 <- reorder(July18_ggheatmap$Gene2, July18_ggheatmap$log2FoldChange)
head(July18_ggheatmap)


ggplot(July18_ggheatmap, aes(x=0, Gene2, fill=log2FoldChange)) +
  geom_tile(aes(width=10), color="black") +
  ylab("Gene ID")+
  xlab("")+
  scale_fill_gradient2(low="cornflowerblue", mid="white", high="firebrick", name="", limits=c(-10,10), labels=c())+
  #guide = guide_legend(label.position="bottom")) +
  coord_flip()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=90, face="bold", size=9), axis.title.x=element_text(size=15, angle=180),
        axis.text.y=element_text(size=0), legend.position = "top", legend.key.width=unit(4,'cm'), legend.text=element_text(size=12, angle=90))
#b=theme(axis.text.x=element_text(angle=90, face="bold", size=9), axis.title.x=element_text(size=15, angle=180),
#        axis.text.y=element_text(size=0), legend.position = "top", legend.key.width=unit(4,'cm'), legend.text=element_text(size=12, angle=90))


#So, what QTL genes are DE at every one of these time points? #make sure these files aren't the whole genome
q <- inner_join(early_v_late_Jun12,early_v_late_Jun29, by="Gene_ID", suffix=c("_early_v_late_Jun12", "_early_v_late_Jun29")) %>%
  dplyr::select("Gene_ID", "log2FoldChange_early_v_late_Jun12", "log2FoldChange_early_v_late_Jun29")
qq <- inner_join(q,zz,by="Gene_ID")
head(qq)
colnames(qq) = c("Gene_ID","June 12", "June 29", "July 18")
head(qq)
#write_csv(qq,"DE_at_every_date_whole_genome.csv")


DE_every_date = qq %>%
  gather(key="Date",value="log2FoldChange",2:4)
DE_every_date$Gene2 <- gsub("-R[A,B]","",DE_every_date$Gene_ID)
DE_every_date$Date <- factor(as.factor(DE_every_date$Date), levels=c("June 12", "June 29", "July 18"))
DE_every_date$Gene2 <- reorder(DE_every_date$Gene2, DE_every_date$log2FoldChange)
head(DE_every_date)

ggplot(DE_every_date, aes(Date, Gene2, fill=log2FoldChange)) +
  geom_tile(aes(), color="black") +
  ylab("Gene ID")+
  xlab("DE Genes at Every Time Point (chr4 QTL region)")+
  scale_fill_gradient2(low="cornflowerblue", mid="white", high="firebrick", name="", limits=c(-10,10), labels=c())+
  coord_flip()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=90, face="bold", size=14), axis.title.x=element_text(size=20, angle=180, face="bold"), axis.title.y=element_text(size=22, face="bold", margin=margin(r=30)),
        axis.text.y=element_text(size=20, angle=90, face="bold", hjust=0.5), legend.position = "top", legend.key.width=unit(4,'cm'), legend.text=element_text(size=12, angle=90))

write_excel_csv(qq,"DEs_every_date_QTL_region.csv")

################

#I'd like to know, in the QTL region:
#1) What genes are differentially expressed at ALL time points that also have genetic variation associated with them?
#incl. modifying variation, high effects, or non-synonymous changes

#What genes are differentially expressed at ANY time point that have genetic variation associated with them?

#cool, now we can figure out what differentially expressed genes also have genetic variants matching the k haplotype, and 
#know what kind of variation it is. #we'll need to start joining lists -- so we gotta reapply the Gene_ID as the first column. 

q <- inner_join(early_v_late_Jun12,early_v_late_Jun29, by="Gene_ID", suffix=c("_early_v_late_Jun12", "_early_v_late_Jun29")) %>%
  dplyr::select("Gene_ID", "log2FoldChange_early_v_late_Jun12", "log2FoldChange_early_v_late_Jun29", "Putative_Function" = "Putative_Function.x_early_v_late_Jun12")
head(q)
zz <- early_v_late_Jul18 %>%
  dplyr::select("Gene_ID", "log2FoldChange_early_v_late_Jul18"="log2FoldChange", "Putative_Function"="Putative_Function.x")
head(zz)

DEGS_every_bloom_date <- inner_join(q,zz, by="Gene_ID") %>%
  dplyr::select("Gene_ID", "June12"="log2FoldChange_early_v_late_Jun12", "June29" = "log2FoldChange_early_v_late_Jun29", "July18"="log2FoldChange_early_v_late_Jul18")
head(DEGS_every_bloom_date)

high_effect <- read.table("num_high_effect_unique_w_counts.txt", header=TRUE)
moderate <- read.table("num_moderate_effect_unique_w_counts.txt", header=TRUE)
modifier_upstream <- read.table("num_UPSTREAM_unique_w_counts.txt", header=TRUE)
modifier_3UTR <- read.table("num_3UTR_unique_w_counts.txt", header=TRUE)
modifier_5UTR <- read.table("num_5UTR_unique_w_counts.txt", header=TRUE)
SVs_high <- read.table("num_transcripts_w_high_SVs_w_counts.txt", header=TRUE)
SVs_upstream <- read.table("num_transcripts_w_upstream_genetic_variants_w_counts.txt", header=TRUE)
sum(str_detect(modifier_3UTR$Gene_ID, '^Pcer_097451')) > 0


##every gene with their different variations
high_effect = high_effect %>%
  dplyr::select("Gene_ID","HIGH_count"="Count")
moderate = moderate %>%
  dplyr::select("Gene_ID","MODERATE_count"="Count")
modifier_3UTR = modifier_3UTR %>%
  dplyr::select("Gene_ID","UTR3_count"="Count")
modifier_5UTR = modifier_5UTR %>%
  dplyr::select("Gene_ID","UTR5_count"="Count")
modifier_upstream = modifier_upstream %>%
  dplyr::select("Gene_ID","UPSTREAM_count"="Count")
high_SV = SVs_high %>%
  dplyr::select("Gene_ID","SVs_HIGH_count"="Count")
upstream_SV = SVs_upstream %>%
  dplyr::select("Gene_ID","SVs_UPSTREAM_count"="Count")

a <- full_join(high_effect,moderate,by="Gene_ID")
b <- full_join(a,modifier_3UTR,by="Gene_ID")
c <- full_join(b,modifier_5UTR,by="Gene_ID")
d <- full_join(c,modifier_upstream,by="Gene_ID")
e <- full_join(d, high_SV, by="Gene_ID")
f <- full_join(e, upstream_SV, by="Gene_ID")
genes_w_variation <- left_join(f,funct_anno,by="Gene_ID")
head(genes_w_variation)
sum(str_detect(genes_w_variation$Gene_ID, '^Pcer_097451')) > 0 #this is my check to know I'm just working with the latest annotation

####### variants for the genes DE at each time point #####
DEGs_w_variation_all_dates = inner_join(genes_w_variation,DEGS_every_bloom_date, by="Gene_ID")
head(DEGs_w_variation_all_dates)
DEGs_w_variation_all_dates= DEGs_w_variation_all_dates %>%
  dplyr::select(Gene_ID, June12,June29,July18,HIGH_count,MODERATE_count,UTR3_count,UTR5_count,UPSTREAM_count,SVs_HIGH_count,SVs_UPSTREAM_count,Putative_Function)
head(DEGs_w_variation_all_dates)
#for venn diagrams
high_effect_DEGS_all = inner_join(high_effect,DEGS_every_bloom_date, by="Gene_ID")
moderate_effect_DEGS_all = inner_join(moderate,DEGS_every_bloom_date, by="Gene_ID")
modifier_5UTR_DEGs_all = inner_join(modifier_5UTR,DEGS_every_bloom_date, by="Gene_ID")
modifier_3UTR_DEGs_all = inner_join(modifier_3UTR,DEGS_every_bloom_date, by="Gene_ID")
modifier_upstream_DEGs_all = inner_join(modifier_upstream,DEGS_every_bloom_date, by="Gene_ID")

#### By date (note these are only for the QTL region)
##June 12th 
head(June12_ggheatmap)
DEGs_w_variation_Jun12 = inner_join(genes_w_variation,June12_ggheatmap, by="Gene_ID")
DEGs_w_variation_Jun12 = DEGs_w_variation_Jun12 %>%
  dplyr::select(Gene_ID,"June12"="log2FoldChange",HIGH_count,MODERATE_count,UTR5_count,UTR3_count,UPSTREAM_count,Putative_Function)
head(DEGs_w_variation_Jun12)
#for venn diagrams
high_effect_DEGS_Jun12 = inner_join(high_effect, June12_ggheatmap,by="Gene_ID")
moderate_effect_DEGS_Jun12 = inner_join(moderate, June12_ggheatmap, by="Gene_ID")
modifier_upstream_DEGs_Jun12 = inner_join(modifier_upstream, June12_ggheatmap, by="Gene_ID")
modifier_5UTR_DEGs_Jun12 = inner_join(modifier_5UTR, June12_ggheatmap, by="Gene_ID")
modifier_3UTR_DEGs_Jun12 = inner_join(modifier_3UTR, June12_ggheatmap, by="Gene_ID")  

##June 29th
head(June29_ggheatmap)
DEGs_w_variation_Jun29 = inner_join(genes_w_variation,June29_ggheatmap, by="Gene_ID")
DEGs_w_variation_Jun29 = DEGs_w_variation_Jun29 %>%
  dplyr::select(Gene_ID,"June29"="log2FoldChange",HIGH_count,MODERATE_count,UTR5_count,UTR3_count,UPSTREAM_count,Putative_Function)
head(DEGs_w_variation_Jun29)
#for venn diagrams
high_effect_DEGS_Jun29 = inner_join(high_effect,June29_ggheatmap, by="Gene_ID")
moderate_effect_DEGS_Jun29 = inner_join(moderate, June29_ggheatmap, by="Gene_ID")
modifier_upstream_DEGs_Jun29 = inner_join(modifier_upstream, June29_ggheatmap, by="Gene_ID")
modifier_5UTR_DEGs_Jun29 = inner_join(modifier_5UTR, June29_ggheatmap, by="Gene_ID")
modifier_3UTR_DEGs_Jun29 = inner_join(modifier_3UTR, June29_ggheatmap, by="Gene_ID")  

#July 18th
head(July18_ggheatmap)
DEGs_w_variation_Jul18 = inner_join(genes_w_variation,July18_ggheatmap, by="Gene_ID")
DEGs_w_variation_Jul18 = DEGs_w_variation_Jul18 %>%
  dplyr::select(Gene_ID,"July18"="log2FoldChange",HIGH_count,MODERATE_count,UTR5_count,UTR3_count,UPSTREAM_count,Putative_Function)
#for venn diagrams
high_effect_DEGS_Jul18 = inner_join(high_effect,July18_ggheatmap, by="Gene_ID")
moderate_effect_DEGS_Jul18 = inner_join(moderate, July18_ggheatmap, by="Gene_ID")
modifier_upstream_DEGs_Jul18 = inner_join(modifier_upstream, July18_ggheatmap, by="Gene_ID")
modifier_5UTR_DEGs_Jul18 = inner_join(modifier_5UTR, July18_ggheatmap, by="Gene_ID")
modifier_3UTR_DEGs_Jul18 = inner_join(modifier_3UTR, July18_ggheatmap, by="Gene_ID") 

DE_Jul18 <- DEGs_w_variation_Jul18 %>%
  dplyr::select(Gene_ID, July18)
DE_June29 <- DEGs_w_variation_Jun29 %>%
  dplyr::select(Gene_ID,June29)
DE_June12 = DEGs_w_variation_Jun12 %>%
  dplyr::select(Gene_ID,June12)

DEGs_any_date_w_variation_a <- full_join(DE_June12,DE_June29, by="Gene_ID")
head(DEGs_any_date_w_variation_a)
DEGs_any_date_w_variation_b <- full_join(DEGs_any_date_w_variation_a,DE_Jul18, by="Gene_ID")
head(DEGs_any_date_w_variation_b)
head(genes_w_variation)
DEGs_any_date_w_variation_c <- left_join(DEGs_any_date_w_variation_b,genes_w_variation,by="Gene_ID")

gene_locations <- read.table("chr4A_genes_w_START.txt", header=TRUE, stringsAsFactors = FALSE) #every chr4A gene with the start coordinate of the mRNA
head(gene_locations)

DEGs_any_date_w_variation_d <- inner_join(DEGs_any_date_w_variation_c,gene_locations,by="Gene_ID") #numbers shouldn't change from c to d if we have all of the genes' starts.

#can also add GO terms and Arabidopsis homologs
homologs <- read.table("sourcherry_to_arabidopsis_hits.txt", header=FALSE, stringsAsFactors = FALSE)
colnames(homologs) = c("Gene_ID", "arabidopsis_blast_hits_homolog")
head(homologs)
homologs <- aggregate(arabidopsis_blast_hits_homolog~Gene_ID,unique(homologs),function(x) paste0(x,collapse = ';'))
head(homologs)
sum(str_detect(homologs$Gene_ID, 'Pcer_097451-RA')) > 0

DEGs_any_date_w_variation_e <- left_join(DEGs_any_date_w_variation_d,homologs,by="Gene_ID")
head(DEGs_any_date_w_variation_e)


####aaaaand GO terms
GO <- read_csv("sourcherry_GOdb_transcripts.csv")
colnames(GO) = c("Gene_ID","GO")
head(GO)

DEGs_any_date_w_variation <- left_join(DEGs_any_date_w_variation_e,GO,by="Gene_ID")
head(DEGs_any_date_w_variation)

write_excel_csv(DEGs_any_date_w_variation,"DEGs_any_date_w_variation.csv") #my precioussss
#any associated of location with type of variation?


#install.packages("VennDiagram")
library(VennDiagram)

venn.diagram(x=list(high_effect_DEGS_all$Gene_ID, moderate_effect_DEGS_all$Gene_ID, 
             modifier_5UTR_DEGs_all$Gene_ID, modifier_3UTR_DEGs_all$Gene_ID, modifier_upstream_DEGs_all$Gene_ID),
             category.names = c("Stop Codon / Frame Shift / Splice Site","AA Change / Codon INDELs", "5'UTR", "3'UTR","Upstream"),
             filename="variant_breakdown_all_dates3.png",
             fill=c("lightpink","violetred3","violetred4","magenta4","lavender"),
             main.fontface="bold",
             sub.fontface="bold",
             lwd=0.4,
             #height = 1000, 
             #width = 1000, 
             resolution = 500,
             compression="lzw",
             margin=0.1,
             main="DEGs at All Dates",
             sub="(with genetic variation by type)",
             main.cex=1.6,
             main.pos=c(0.5, 1.05),
             sub.cex=1.2,
             sub.pos=c(0.5, 1),
             cat.dist=c(0.18, 0.20, 0.23, 0.22, 0.205),
             cat.fontface = "bold",
             cat.default.pos="outer",
             cat.pos=c(0,345,270,90,350),
             rotation.degree=355,
             cat.cex=0.88)

#there must be some kind of limit... like 5 entries or something? Guess I'll just look at these 5, variants landing in genes. 

#number of genes differentially expressed, and what each date has in common. (again, QTL region)
venn.diagram(x=list(early_v_late_Jun12$Gene_ID, early_v_late_Jun29$Gene_ID, early_v_late_Jul18$Gene_ID),
             category.names = c("June 12","June 29", "July 18"),
             filename="shared_DEGS_by_date_in_QTL_region.png",
             fill=c("forestgreen","gold","skyblue"),
             main.fontface="bold",
             sub.fontface="bold",
             lwd=0.4,
             #height = 1000, 
             #width = 1000, 
             resolution = 700,
             compression="lzw",
             margin=0.04,
             cat.fontface = "bold",
             cat.fontsize = 10,
             main="Number of DEGs by Date",
             main.pos=c(0.5, 1.05),
             sub="(chr4 QTL region)",
             sub.pos=c(0.5, 1),
             main.cex=1.4,
             sub.cex=1.0,
             )

###NOTE- you must overwrite the current early_v_late_date files above before getting numbers for the full genome:
#number of genes differentially expressed, and what each date has in common. (again, QTL region)
venn.diagram(x=list(early_v_late_Jun12$Gene_ID, early_v_late_Jun29$Gene_ID, early_v_late_Jul18$Gene_ID),
             category.names = c("June 12","June 29", "July 18"),
             filename="shared_DEGS_by_date_whole_genome.png",
             fill=c("forestgreen","gold","skyblue"),
             main.fontface="bold",
             sub.fontface="bold",
             lwd=0.4,
             #ext.line.lty=0.01,
             #height = 1000, 
             #width = 1000, 
             resolution = 700,
             compression="lzw",
             margin=0.04,
             cat.fontface = "bold",
             cat.fontsize = 10,
             main="Number of DEGs by Date",
             main.pos=c(0.5, 1.05),
             sub="(full genome)",
             sub.pos=c(0.5, 1),
             main.cex=1.4,
             sub.cex=1.0,
)


venn.diagram(x=list(high_effect_DEGS_Jun12$Gene_ID, moderate_effect_DEGS_Jun12$Gene_ID, 
                    modifier_5UTR_DEGs_Jun12$Gene_ID, modifier_3UTR_DEGs_Jun12$Gene_ID, modifier_upstream_DEGs_Jun12$Gene_ID),
             category.names = c("Stop Codon / Frame Shift / Splice Site","AA Change / Codon INDELs", "5'UTR", "3'UTR","Upstream"),
             filename="variant_breakdown_June12.png",
             fill=c("darkgreen","green4","chartreuse3","darkolivegreen2","darkolivegreen1"),
             main.fontface="bold",
             sub.fontface="bold",
             lwd=0.4,
             #height = 1000, 
             #width = 1000, 
             resolution = 500,
             compression="lzw",
             margin=0.1,
             main="DEGs on June 12th",
             sub="(with genetic variation by type)",
             main.cex=1.6,
             main.pos=c(0.5, 1.05),
             sub.cex=1.2,
             sub.pos=c(0.5, 1),
             cat.dist=c(0.18, 0.20, 0.23, 0.22, 0.205),
             cat.fontface = "bold",
             cat.default.pos="outer",
             cat.pos=c(0,345,270,90,350),
             rotation.degree=355,
             cat.cex=0.88)

venn.diagram(x=list(high_effect_DEGS_Jun29$Gene_ID, moderate_effect_DEGS_Jun29$Gene_ID, 
                    modifier_5UTR_DEGs_Jun29$Gene_ID, modifier_3UTR_DEGs_Jun29$Gene_ID, modifier_upstream_DEGs_Jun29$Gene_ID),
             category.names = c("Stop Codon / Frame Shift / Splice Site","AA Change / Codon INDELs", "5'UTR", "3'UTR","Upstream"),
             filename="variant_breakdown_Jun29.png",
             fill=c("orange2","orange1","goldenrod2","gold2","gold"),
             main.fontface="bold",
             sub.fontface="bold",
             lwd=0.4,
             #height = 1000, 
             #width = 1000, 
             resolution = 500,
             compression="lzw",
             margin=0.1,
             main="DEGs on June 29th",
             sub="(with genetic variation by type)",
             main.cex=1.6,
             main.pos=c(0.5, 1.05),
             sub.cex=1.2,
             sub.pos=c(0.5, 1),
             cat.dist=c(0.18, 0.20, 0.23, 0.22, 0.205),
             cat.fontface = "bold",
             cat.default.pos="outer",
             cat.pos=c(0,345,270,90,350),
             rotation.degree=355,
             cat.cex=0.88)

venn.diagram(x=list(high_effect_DEGS_Jul18$Gene_ID, moderate_effect_DEGS_Jul18$Gene_ID, 
                    modifier_5UTR_DEGs_Jul18$Gene_ID, modifier_3UTR_DEGs_Jul18$Gene_ID, modifier_upstream_DEGs_Jul18$Gene_ID),
             category.names = c("Stop Codon / Frame Shift / Splice Site","AA Change / Codon INDELs", "5'UTR", "3'UTR","Upstream"),
             filename="variant_breakdown_Jul18.png",
             fill=c("deepskyblue4","skyblue3","skyblue","lightskyblue1","lightblue1"),
             main.fontface="bold",
             sub.fontface="bold",
             lwd=0.4,
             #height = 1000, 
             #width = 1000, 
             resolution = 500,
             compression="lzw",
             margin=0.1,
             main="DEGs on July 18th",
             sub="(with genetic variation by type)",
             main.cex=1.6,
             main.pos=c(0.5, 1.05),
             sub.cex=1.2,
             sub.pos=c(0.5, 1),
             cat.dist=c(0.18, 0.20, 0.23, 0.22, 0.205),
             cat.fontface = "bold",
             cat.default.pos="outer",
             cat.pos=c(0,345,270,90,350),
             rotation.degree=355,
             cat.cex=0.88)

#DEGs and select genetic variation
head(high_effect_DEGS_all)

#all_dates
write_excel_csv(high_effect_DEGS_all,"high_effect_DEGs_QTL_region_every_date.csv")
write_excel_csv(moderate_effect_DEGS_all, "nonsynonymous_effect_DEGs_QTL_region_every_date.csv")
write_excel_csv(modifier_5UTR_DEGs_all, "5UTR_DEGs_QTL_region_every_date.csv")
write_excel_csv(modifier_3UTR_DEGs_all, "3UTR_DEGs_QTL_region_every_date.csv")
write_excel_csv(modifier_upstream_DEGs_all, "upstream_DEGs_QTL_region_every_date.csv")
#Jun12
write_excel_csv(high_effect_DEGS_Jun12, "high_effect_DEGs_QTL_region_Jun12.csv")
write_excel_csv(moderate_effect_DEGS_Jun12, "nonsyn_DEGs_QTL_region_Jun12.csv")
write_excel_csv(modifier_upstream, "upstream_modifier_DEGs_QTL_region_Jun12.csv")
write_excel_csv(modifier_5UTR_DEGs_Jun12, "5UTR_modifier_DEGs_QTL_region_Jun12.csv")
write_excel_csv(modifier_3UTR_DEGs_Jun12, "3UTR_modifier_DEGs_QTL_region_Jun12.csv")
#Jun29
write_excel_csv(high_effect_DEGS_Jun29, "high_effect_DEGs_QTL_region_Jun29.csv")
write_excel_csv(moderate_effect_DEGS_Jun29, "nonsyn_DEGs_QTL_region_Jun29.csv")
write_excel_csv(modifier_upstream, "upstream_modifier_DEGs_QTL_region_Jun29.csv")
write_excel_csv(modifier_5UTR_DEGs_Jun29, "5UTR_modifier_DEGs_QTL_region_Jun29.csv")
write_excel_csv(modifier_3UTR_DEGs_Jun29, "3UTR_modifier_DEGs_QTL_region_Jun29.csv")
#Jul18
write_excel_csv(high_effect_DEGS_Jul18, "high_effect_DEGs_QTL_region_Jul18.csv")
write_excel_csv(moderate_effect_DEGS_Jul18, "nonsyn_DEGs_QTL_region_Jul18.csv")
write_excel_csv(modifier_upstream, "upstream_modifier_DEGs_QTL_region_Jul18.csv")
write_excel_csv(modifier_5UTR_DEGs_Jul18, "5UTR_modifier_DEGs_QTL_region_Jul18.csv")
write_excel_csv(modifier_3UTR_DEGs_Jul18, "3UTR_modifier_DEGs_QTL_region_Jul18.csv")




##################################################################################################################
#                                                                                                                #  
# GO enrichment for the differentially expressed genes of early- and late-bloomers around the time of initiation #
#                                                                                                                #
##################################################################################################################

#see the end of the WGCNA script to see how our GO database was created. It's 
#a combo of interproscan GO terms for sour cherry
#and GO terms for the reciprocal blast hits of arabidopsis and sour cherry.
#See scripts for blast in this folder for the parameters used, and the 'steps 2 get GO terms w arabidopsis . txt" to 
#view some commands that were used to create some of the 'pre-files' of the database.

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/DESeq2/QTL_round2_update/")

library(tidyverse)
library(dplyr)
library(topGO)
#BiocManager::install("apeglm")
library(apeglm)
source("topGOFunctions.R") #Ben kindly lent us his function wrapper that he made for his cucumber paper.

#prepare as named list for use with topGO
GOdatabase <- read.csv("sourcherry_GOdb_transcripts.csv", header=TRUE, stringsAsFactors = FALSE)
#create a list where GO terms are named by their gene name
GOList<-setNames(nm = GOdatabase$GENE, strsplit(GOdatabase$GO, ";"))
head(GOList)

list2env(dds.res.flt2,envir=.GlobalEnv) #NOTE - this overwrites the variables to be the results of the whole genome, not just
#the QTL region. 


#remember that late is the reference - so genes that are up are higher in the earlies, and those that are down are lower in the earlies 
early_v_late_Jun12_UP <- early_v_late_Jun12[early_v_late_Jun12$log2FoldChange > 0, ]
early_v_late_Jun12_DOWN <- early_v_late_Jun12[early_v_late_Jun12$log2FoldChange < 0, ]
early_v_late_Jun29_UP <- early_v_late_Jun29[early_v_late_Jun29$log2FoldChange > 0, ]
early_v_late_Jun29_DOWN <- early_v_late_Jun29[early_v_late_Jun29$log2FoldChange < 0, ]
early_v_late_Jul18_UP <- early_v_late_Jul18[early_v_late_Jul18$log2FoldChange > 0, ]
early_v_late_Jul18_DOWN <- early_v_late_Jul18[early_v_late_Jul18$log2FoldChange < 0, ]
early_Jun12_Jun29_UP <- early_Jun12_Jun29[early_Jun12_Jun29$log2FoldChange > 0, ]
early_Jun12_Jun29_DOWN <- early_Jun12_Jun29[early_Jun12_Jun29$log2FoldChange < 0, ]
early_Jun29_Jul18_UP <- early_Jun29_Jul18[early_Jun29_Jul18$log2FoldChange > 0, ]
early_Jun29_Jul18_DOWN <- early_Jun29_Jul18[early_Jun29_Jul18$log2FoldChange < 0, ]
late_Jun12_Jun29_UP <- late_Jun12_Jun29[late_Jun12_Jun29$log2FoldChange > 0, ]
late_Jun12_Jun29_DOWN <- late_Jun12_Jun29[late_Jun12_Jun29$log2FoldChange < 0, ]
late_Jun29_Jul18_UP <- late_Jun29_Jul18[late_Jun29_Jul18$log2FoldChange > 0, ]
late_Jun29_Jul18_DOWN <- late_Jun29_Jul18[late_Jun29_Jul18$log2FoldChange < 0, ]


#I'd like to know what sort of GO terms are enriched in each of the above categories
#BP = Biological Process
#info from the manual:
#The elim and weight algorithms were introduced in Alexa et al. (2006). The default algorithm used by the
#topGO package is a mixture between the elim and the weight algorithms and it will be referred as weight01.

#####Jun12
#UP
GO_early_v_late_Jun12_UP <- runTopGoAnalysis(DEgeneSet = early_v_late_Jun12_UP$Gene_ID,
                            dds = subA_dds,
                            GOdb = GOList,
                            onts = "BP",
                            nodeSize = 50)

top50_enriched_early_v_late_Jun12_UP <- exportGOtable(GO_early_v_late_Jun12_UP, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))%>%
  mutate(Direction="UP")
head(top50_enriched_early_v_late_Jun12_UP)

#DOWN
GO_early_v_late_Jun12_DOWN <- runTopGoAnalysis(DEgeneSet = early_v_late_Jun12_DOWN$Gene_ID,
                                             dds = subA_dds,
                                             GOdb = GOList,
                                             onts = "BP",
                                             nodeSize = 50)
#BP = Biological Process
top50_enriched_early_v_late_Jun12_DOWN <- exportGOtable(GO_early_v_late_Jun12_DOWN, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01)) %>%
  mutate(Direction="DOWN")
head(top50_enriched_early_v_late_Jun12_DOWN)
top50_enriched_GOterms_early_v_late_Jun12 <- rbind(top50_enriched_early_v_late_Jun12_UP,top50_enriched_early_v_late_Jun12_DOWN)


#####Jun29
#UP
GO_early_v_late_Jun29_UP <- runTopGoAnalysis(DEgeneSet = early_v_late_Jun29_UP$Gene_ID,
                                             dds = subA_dds,
                                             GOdb = GOList,
                                             onts = "BP",
                                             nodeSize = 50)


top50_enriched_early_v_late_Jun29_UP <- exportGOtable(GO_early_v_late_Jun29_UP, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))%>%
  mutate(Direction="UP")
head(top50_enriched_early_v_late_Jun29_UP)

#DOWN
GO_early_v_late_Jun29_DOWN <- runTopGoAnalysis(DEgeneSet = early_v_late_Jun29_DOWN$Gene_ID,
                                             dds = subA_dds,
                                             GOdb = GOList,
                                             onts = "BP",
                                             nodeSize = 50)


top50_enriched_early_v_late_Jun29_DOWN <- exportGOtable(GO_early_v_late_Jun29_DOWN, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))%>%
  mutate(Direction="DOWN")
head(top50_enriched_early_v_late_Jun29_DOWN)

top50_enriched_GOterms_early_v_late_Jun29 <- rbind(top50_enriched_early_v_late_Jun29_UP,top50_enriched_early_v_late_Jun29_DOWN)


#####Jul18
#UP
GO_early_v_late_Jul18_UP <- runTopGoAnalysis(DEgeneSet = early_v_late_Jul18_UP$Gene_ID,
                                             dds = subA_dds,
                                             GOdb = GOList,
                                             onts = "BP",
                                             nodeSize = 50)


top50_enriched_early_v_late_Jul18_UP <- exportGOtable(GO_early_v_late_Jul18_UP, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01)) %>%
  mutate(Direction = "UP")
head(top50_enriched_early_v_late_Jul18_UP)

#DOWN
GO_early_v_late_Jul18_DOWN <- runTopGoAnalysis(DEgeneSet = early_v_late_Jul18_DOWN$Gene_ID,
                                               dds = subA_dds,
                                               GOdb = GOList,
                                               onts = "BP",
                                               nodeSize = 50)


top50_enriched_early_v_late_Jul18_DOWN <- exportGOtable(GO_early_v_late_Jul18_DOWN, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))%>%
  mutate(Direction = "DOWN")

head(top50_enriched_early_v_late_Jul18_DOWN)


top50_enriched_GOterms_early_v_late_Jul18 <- rbind(top50_enriched_early_v_late_Jul18_UP,top50_enriched_early_v_late_Jul18_DOWN)

#Positive regulation of flower development is an enriched GO term for DEGs more highly expressed in the early-bloomers on Jul18
#On the other hand, gibberellin biosynthetic process is 
#That's reassuring, lol.

#I don't care as much about these comparisons, but a reviewer might ask.

#early from Jun12 to Jun29
#UP - means it's higher on Jun29 and these genes increased in expression from Jun12 to Jun29
GO_early_Jun12_Jun29_UP <- runTopGoAnalysis(DEgeneSet = early_Jun12_Jun29_UP$Gene_ID,
                                            dds = subA_dds,
                                            GOdb = GOList,
                                            onts = "BP",
                                            nodeSize = 50)


top50_enriched_early_Jun12_Jun29_UP <- exportGOtable(GO_early_Jun12_Jun29_UP, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01)) %>%
  mutate(Direction="UP")

head(top50_enriched_early_Jun12_Jun29_UP)

#DOWN - means it was higher on June12 and these genes decreased in expression by Jun29.
GO_early_Jun12_Jun29_DOWN <- runTopGoAnalysis(DEgeneSet = early_Jun12_Jun29_DOWN$Gene_ID,
                                              dds = subA_dds,
                                              GOdb = GOList,
                                              onts = "BP",
                                              nodeSize = 50)


top50_enriched_early_Jun12_Jun29_DOWN <- exportGOtable(GO_early_Jun12_Jun29_DOWN, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))%>%
  mutate(Direction="DOWN")

head(top50_enriched_early_Jun12_Jun29_DOWN)

top50_enriched_GOterms_early_Jun12_Jun29 <- rbind(top50_enriched_early_Jun12_Jun29_UP,top50_enriched_early_Jun12_Jun29_DOWN)


#early from Jun29 to Jul18
#UP - means it's higher on July18 and these genes increased in expression from Jun29 to July18
GO_early_Jun29_Jul18_UP <- runTopGoAnalysis(DEgeneSet = early_Jun29_Jul18_UP$Gene_ID,
                                            dds = subA_dds,
                                            GOdb = GOList,
                                            onts = "BP",
                                            nodeSize = 50)


top50_enriched_early_Jun29_Jul18_UP <- exportGOtable(GO_early_Jun29_Jul18_UP, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01)) %>%
  mutate(Direction = "UP")

head(top50_enriched_early_Jun29_Jul18_UP)

#DOWN - means it was higher on June29 and these genes decreased in expression by July18.
GO_early_Jun29_Jul18_DOWN <- runTopGoAnalysis(DEgeneSet = early_Jun29_Jul18_DOWN$Gene_ID,
                                              dds = subA_dds,
                                              GOdb = GOList,
                                              onts = "BP",
                                              nodeSize = 50)


top50_enriched_early_Jun29_Jul18_DOWN <- exportGOtable(GO_early_Jun29_Jul18_DOWN, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01)) %>%
  mutate(Direction = "DOWN")

head(top50_enriched_early_Jun29_Jul18_DOWN)

top50_enriched_GOterms_early_Jun29_Jul18 <- rbind(top50_enriched_early_Jun29_Jul18_UP,top50_enriched_early_Jun29_Jul18_DOWN)



################Now let's look at the late-bloomers' journey

#late from Jun12 to Jun29
#UP - means it's higher on Jun29 and these genes increased in expression from Jun12 to Jun29
GO_late_Jun12_Jun29_UP <- runTopGoAnalysis(DEgeneSet = late_Jun12_Jun29_UP$Gene_ID,
                                            dds = subA_dds,
                                            GOdb = GOList,
                                            onts = "BP",
                                            nodeSize = 50)


top50_enriched_late_Jun12_Jun29_UP <- exportGOtable(GO_late_Jun12_Jun29_UP, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))%>%
  mutate(Direction="UP")

#DOWN - means it was higher on June12 and these genes decreased in expression by Jun29.
GO_late_Jun12_Jun29_DOWN <- runTopGoAnalysis(DEgeneSet = late_Jun12_Jun29_DOWN$Gene_ID,
                                              dds = subA_dds,
                                              GOdb = GOList,
                                              onts = "BP",
                                              nodeSize = 50)


top50_enriched_late_Jun12_Jun29_DOWN <- exportGOtable(GO_late_Jun12_Jun29_DOWN, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))%>%
  mutate(Direction="DOWN")

top50_enriched_GOterms_late_Jun12_Jun29 <- rbind(top50_enriched_late_Jun12_Jun29_UP,top50_enriched_late_Jun12_Jun29_DOWN)


#late from Jun29 to Jul18
#UP - means it's higher on July18 and these genes increased in expression from Jun29 to July18
GO_late_Jun29_Jul18_UP <- runTopGoAnalysis(DEgeneSet = late_Jun29_Jul18_UP$Gene_ID,
                                            dds = subA_dds,
                                            GOdb = GOList,
                                            onts = "BP",
                                            nodeSize = 50)


top50_enriched_late_Jun29_Jul18_UP <- exportGOtable(GO_late_Jun29_Jul18_UP, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01)) %>%
  mutate(Direction = "UP")

#DOWN - means it was higher on June29 and these genes decreased in expression by July18.
GO_late_Jun29_Jul18_DOWN <- runTopGoAnalysis(DEgeneSet = late_Jun29_Jul18_DOWN$Gene_ID,
                                              dds = subA_dds,
                                              GOdb = GOList,
                                              onts = "BP",
                                              nodeSize = 50)


top50_enriched_late_Jun29_Jul18_DOWN <- exportGOtable(GO_late_Jun29_Jul18_DOWN, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01)) %>%
  mutate(Direction = "DOWN")

top50_enriched_GOterms_late_Jun29_Jul18 <- rbind(top50_enriched_late_Jun29_Jul18_UP,top50_enriched_late_Jun29_Jul18_DOWN)

###############what about the differentially expressed genes that are always up or always down at all dates?

UP_at_all_dates <- inner_join(early_v_late_Jun12_UP, early_v_late_Jun29_UP, by="Gene_ID") %>%
  inner_join(., early_v_late_Jul18_UP, by="Gene_ID")

DOWN_at_all_dates <- inner_join(early_v_late_Jun12_DOWN, early_v_late_Jun29_DOWN, by="Gene_ID") %>%
  inner_join(., early_v_late_Jul18_DOWN, by="Gene_ID")

###########topGO

GO_UP_at_all_dates <- runTopGoAnalysis(DEgeneSet = UP_at_all_dates$Gene_ID,
                                         dds = subA_dds,
                                         GOdb = GOList,
                                         onts = "BP",
                                         nodeSize = 50)


top50_enriched_UP_at_all_dates <- exportGOtable(GO_UP_at_all_dates, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01)) %>%
  mutate(Direction = "UP")

head(top50_enriched_UP_at_all_dates)

#########

GO_DOWN_at_all_dates <- runTopGoAnalysis(DEgeneSet = DOWN_at_all_dates$Gene_ID,
                                           dds = subA_dds,
                                           GOdb = GOList,
                                           onts = "BP",
                                           nodeSize = 50)


top50_enriched_DOWN_at_all_dates <- exportGOtable(GO_DOWN_at_all_dates, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01)) %>%
  mutate(Direction = "DOWN")

head(top50_enriched_DOWN_at_all_dates)

top50_enriched_DEGs_at_all_dates <- rbind(top50_enriched_UP_at_all_dates, top50_enriched_DOWN_at_all_dates)



write_tsv(top50_enriched_GOterms_early_v_late_Jun12, "top50_enriched_GOterms_early_v_late_Jun12.tsv")
write_tsv(top50_enriched_GOterms_early_v_late_Jun29, "top50_enriched_GOterms_early_v_late_Jun29.tsv")
write_tsv(top50_enriched_GOterms_early_v_late_Jul18, "top50_enriched_GOterms_early_v_late_Jul18.tsv")
write_tsv(top50_enriched_GOterms_early_Jun12_Jun29, "top50_enriched_GOterms_early_Jun12_Jun29.tsv")
write_tsv(top50_enriched_GOterms_early_Jun29_Jul18, "top50_enriched_GOterms_early_Jun29_Jul18.tsv" )
write_tsv(top50_enriched_GOterms_late_Jun12_Jun29, "top50_enriched_GOterms_late_Jun12_Jun29.tsv")
write_tsv(top50_enriched_GOterms_late_Jun29_Jul18, "top50_enriched_GOterms_late_Jun29_Jul18.tsv")
write_tsv(top50_enriched_DEGs_at_all_dates, "top50_enriched_DEGs_at_all_dates.tsv")

######################lastly, what about the genes specifically in the hotpink module that we discovered in our WGCNA analysis?
######And other modules too!

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/WGCNA/post Manfeld:Goeckeritz meeting/gene_atlas_flt5_cov0.2_modules_cut_default/")

#note that, on average, this module is defined as most of these genes having higher expression in the early-bloomers. 

#prepare as named list for use with topGO -- note we need to change it to gene name, not transcript name. 
#GOdatabase <- read.csv("sourcherry_GOdb_transcripts.csv", header=TRUE, stringsAsFactors = FALSE)
GOdatabase$GENE<- gsub("-R[A,B]","",fixed = FALSE, GOdatabase$GENE)
head(GOdatabase)

#create a list where GO terms are named by their gene name
GOList<-setNames(nm = GOdatabase$GENE, strsplit(GOdatabase$GO, ";"))
head(GOList)

##I'm too lazy to create a loop XD
hotpink <- read.table("hotpink_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
violetred3 <- read.table("violetred3_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
black <- read.table("black_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
chartreuse4 <- read.table("chartreuse4_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
chocolate2 <- read.table("chocolate2_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
darkgreen <- read.table("darkgreen_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
dodgerblue3 <- read.table("dodgerblue3_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
firebrick3 <- read.table("firebrick3_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
gold <- read.table("gold_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
khaki <- read.table("khaki_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
seashell <- read.table("seashell_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
lavender <- read.table("lavender_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
lightblue2 <- read.table("lightblue2_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
lightcoral <- read.table("lightcoral_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
pink <- read.table("pink_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
palevioletred1 <- read.table("palevioletred1_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
plum4 <- read.table("plum4_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
yellowgreen <- read.table("yellowgreen_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
grey <- read.table("grey_module_genes.txt", header=FALSE, stringsAsFactors=FALSE)
royalblue4 <- read.table("royalblue4_module_genes.txt", header=FALSE, stringsAsFactors = FALSE)

#oh shoot, I don't think I can use Ben's script here -- unless I do something about the gene vs transcript names in the dds object. 
#so, I could just run the DESeq commands again to get the dds object to contain gene names, not transcript names. 
#the results will be virtually the same since only a handful of genes (< 20) contain secondary transcripts (-RB).


################################ Prep a dds object with gene names, not transcript names, for use in WGCNA.


setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/DESeq2/QTL_round2_update/")
subA_gene_countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))

colnames(subA_gene_countData) <- sub("^X","",perl = T, colnames(subA_gene_countData))
colnames(subA_gene_countData) <- gsub(".","_",fixed = TRUE, colnames(subA_gene_countData))
head(subA_gene_countData)
subA_countsNonZero_gene <- subA_gene_countData[apply(subA_gene_countData,1,function(x){!all(x == 0)}),]
nrow(subA_countsNonZero_gene)
subA_countsNonZero_gene_df <- as.data.frame(subA_countsNonZero_gene)
subA_countsNonZero_gene_no_outlier <- as.matrix(subA_countsNonZero_gene_df[,!names(subA_countsNonZero_gene_df) %in% c("27_04_34_6_12_19_R2","27_03_46_6_29_19_R3")])
head(subA_countsNonZero_gene_no_outlier)
subA_coldata_gene_no_outlier <- subA_coldata[!row.names(subA_coldata) %in% c("27_04_34_6_12_19_R2","27_03_46_6_29_19_R3"),]
all(rownames(subA_coldata_gene_no_outlier) == colnames(subA_countsNonZero_gene_no_outlier))
empty_rows = data.frame(matrix(ncol=ncol(subA_countsNonZero_gene_no_outlier), nrow=3))
colnames(empty_rows) = colnames(subA_countsNonZero_gene_no_outlier)
subA_w_factors = rbind(subA_countsNonZero_gene_no_outlier, empty_rows)

for (i in 1:ncol(subA_w_factors)){
  if (grepl("27_02_19",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_02_19"
  } else if (grepl("27_03_08",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_03_08"
  } else if (grepl("27_03_25",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_03_25"
  } else if (grepl("27_03_46",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_03_46"
  } else if (grepl("27_04_12",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_04_12"
  } else if (grepl("27_04_34",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)-2),i] <- "27_04_34"
  }
  if (grepl("6_12_19",colnames(subA_w_factors)[i])){
    subA_w_factors[nrow(subA_w_factors)-1,i] <- "Jun12"
  } else if (grepl("6_29_19",colnames(subA_w_factors)[i])){
    subA_w_factors[nrow(subA_w_factors)-1,i] <- "Jun29"
  }  else if (grepl("7_18_19",colnames(subA_w_factors)[i])){
    subA_w_factors[nrow(subA_w_factors)-1,i] <- "Jul18"
  }
  if (grepl("27_02_19",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "late"
  } else if (grepl("27_03_08",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "late"
  } else if (grepl("27_03_25",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "early"
  } else if (grepl("27_03_46",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "early"
  } else if (grepl("27_04_12",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "late"
  } else if (grepl("27_04_34",colnames(subA_w_factors)[i])){
    subA_w_factors[(nrow(subA_w_factors)),i] <- "early"
  }
}

#subA_coldata was made earlier in the script wo we're good. 

#to prep individual tree to be nested in our linear model. 
for (i in 1:nrow(subA_coldata_gene_no_outlier)){
  if (grepl("27_02_19", subA_coldata_gene_no_outlier$tree[i])){
    subA_coldata_gene_no_outlier[i, "nested_individual"] <- "1"
  } else if (grepl("27_03_08", subA_coldata_gene_no_outlier$tree[i])){
    subA_coldata_gene_no_outlier[i, "nested_individual"] <- "2"
  } else if (grepl("27_04_12", subA_coldata_gene_no_outlier$tree[i])){
    subA_coldata_gene_no_outlier[i, "nested_individual"] <- "3"
  } else if (grepl("27_03_25", subA_coldata_gene_no_outlier$tree[i])){
    subA_coldata_gene_no_outlier[i, "nested_individual"] <- "1"
  } else if (grepl("27_03_46", subA_coldata_gene_no_outlier$tree[i])){
    subA_coldata_gene_no_outlier[i, "nested_individual"] <- "2"
  } else if (grepl("27_04_34", subA_coldata_gene_no_outlier$tree[i])){
    subA_coldata_gene_no_outlier[i, "nested_individual"] <- "3"
  }
}

all(rownames(subA_coldata_gene_no_outlier) == colnames(subA_countsNonZero_gene_no_outlier))

subA_dds1 <- DESeqDataSetFromMatrix(subA_countsNonZero_gene_no_outlier, 
                                     colData=subA_coldata_gene_no_outlier, 
                                     design= ~ bloom + bloom:nested_individual + bloom:date + date)

subA_dds2 = estimateSizeFactors(subA_dds1)
subA_dds2 = estimateDispersions(subA_dds2)
subA_dds2 = nbinomWaldTest(subA_dds2, maxit=100000)

resultsNames(subA_dds2) 

subA_dds2$bloom <- relevel(subA_dds2$bloom,ref = "late") 
subA_dds2$date <- relevel(subA_dds2$date,ref="Jun12") 
subA_dds2$group <- factor(paste0(subA_dds2$bloom,subA_dds2$date)) 
design(subA_dds2) <- ~group
head(subA_dds2)

#since we regrouped the variables, we need to estimate the variance again:
subA_dds = DESeq(subA_dds2) ###DESeq object for WGCNA's GO enrichment
str(subA_dds2)
subA_dds2$group


########### coolio, we in business boizzzz

#send the output files to the WGCNA folder
setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/WGCNA/post Manfeld:Goeckeritz meeting/")


#hotpink
GO_hotpink <- runTopGoAnalysis(DEgeneSet = hotpink$V1,
                                         dds = subA_dds2,
                                         GOdb = GOList,
                                         onts = "BP",
                                         nodeSize = 50)


top50_enriched_hotpink <- exportGOtable(GO_hotpink, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

head(top50_enriched_hotpink)

write_tsv(top50_enriched_hotpink, "top50_enriched_hotpink.tsv")

#violetred3
GO_violetred3 <- runTopGoAnalysis(DEgeneSet = violetred3$V1,
                                  dds = subA_dds2,
                                  GOdb = GOList,
                                  onts = "BP",
                                  nodeSize = 50)


top50_enriched_violetred3 <- exportGOtable(GO_violetred3, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

head(top50_enriched_violetred3)

write_tsv(top50_enriched_violetred3, "top50_enriched_violetred3.tsv")

#black
GO_black <- runTopGoAnalysis(DEgeneSet = black$V1,
                                  dds = subA_dds2,
                                  GOdb = GOList,
                                  onts = "BP",
                                  nodeSize = 50)


top50_enriched_black <- exportGOtable(GO_black, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

head(top50_enriched_black)

write_tsv(top50_enriched_black, "top50_enriched_black.tsv")

#chartreuse4
GO_chartreuse4 <- runTopGoAnalysis(DEgeneSet = chartreuse4$V1,
                             dds = subA_dds2,
                             GOdb = GOList,
                             onts = "BP",
                             nodeSize = 50)


top50_enriched_chartreuse4 <- exportGOtable(GO_chartreuse4, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

head(top50_enriched_chartreuse4)

write_tsv(top50_enriched_chartreuse4, "top50_enriched_chartreuse4.tsv")

#chocolate2
GO_chocolate2 <- runTopGoAnalysis(DEgeneSet = chocolate2$V1,
                                   dds = subA_dds2,
                                   GOdb = GOList,
                                   onts = "BP",
                                   nodeSize = 50)


top50_enriched_chocolate2 <- exportGOtable(GO_chocolate2, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

head(top50_enriched_chocolate2)

write_tsv(top50_enriched_chocolate2, "top50_enriched_chocolate2.tsv")

#darkgreen
GO_darkgreen <- runTopGoAnalysis(DEgeneSet = darkgreen$V1,
                                  dds = subA_dds2,
                                  GOdb = GOList,
                                  onts = "BP",
                                  nodeSize = 50)


top50_enriched_darkgreen <- exportGOtable(GO_darkgreen, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

head(top50_enriched_darkgreen)

write_tsv(top50_enriched_darkgreen, "top50_enriched_darkgreen.tsv")

#dodgerblue3
GO_dodgerblue3 <- runTopGoAnalysis(DEgeneSet = dodgerblue3$V1,
                                 dds = subA_dds2,
                                 GOdb = GOList,
                                 onts = "BP",
                                 nodeSize = 50)


top50_enriched_dodgerblue3 <- exportGOtable(GO_dodgerblue3, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

head(top50_enriched_dodgerblue3)

write_tsv(top50_enriched_dodgerblue3, "top50_enriched_dodgerblue3.tsv")

#firebrick3
GO_firebrick3 <- runTopGoAnalysis(DEgeneSet = firebrick3$V1,
                                   dds = subA_dds2,
                                   GOdb = GOList,
                                   onts = "BP",
                                   nodeSize = 50)


top50_enriched_firebrick3 <- exportGOtable(GO_firebrick3, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_firebrick3, "top50_enriched_firebrick3.tsv")

head(top50_enriched_firebrick3)

#gold
GO_gold <- runTopGoAnalysis(DEgeneSet = gold$V1,
                                  dds = subA_dds2,
                                  GOdb = GOList,
                                  onts = "BP",
                                  nodeSize = 50)


top50_enriched_gold <- exportGOtable(GO_gold, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_gold, "top50_enriched_gold.tsv")

#lavender
GO_lavender <- runTopGoAnalysis(DEgeneSet = lavender$V1,
                            dds = subA_dds2,
                            GOdb = GOList,
                            onts = "BP",
                            nodeSize = 50)


top50_enriched_lavender <- exportGOtable(GO_lavender, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_lavender, "top50_enriched_lavender.tsv")

#khaki
GO_khaki <- runTopGoAnalysis(DEgeneSet = khaki$V1,
                                dds = subA_dds2,
                                GOdb = GOList,
                                onts = "BP",
                                nodeSize = 50)


top50_enriched_khaki <- exportGOtable(GO_khaki, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_khaki, "top50_enriched_khaki.tsv")

#lightblue2
GO_lightblue2 <- runTopGoAnalysis(DEgeneSet = lightblue2$V1,
                             dds = subA_dds2,
                             GOdb = GOList,
                             onts = "BP",
                             nodeSize = 50)


top50_enriched_lightblue2 <- exportGOtable(GO_lightblue2, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_lightblue2, "top50_enriched_lightblue2.tsv")

#lightcoral
GO_lightcoral <- runTopGoAnalysis(DEgeneSet = lightcoral$V1,
                                  dds = subA_dds2,
                                  GOdb = GOList,
                                  onts = "BP",
                                  nodeSize = 50)


top50_enriched_lightcoral <- exportGOtable(GO_lightcoral, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_lightcoral, "top50_enriched_lightcoral.tsv")

#palevioletred1
GO_palevioletred1 <- runTopGoAnalysis(DEgeneSet = palevioletred1$V1,
                                  dds = subA_dds2,
                                  GOdb = GOList,
                                  onts = "BP",
                                  nodeSize = 50)


top50_enriched_palevioletred1 <- exportGOtable(GO_palevioletred1, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_palevioletred1, "top50_enriched_palevioletred1.tsv")

#seashell
GO_seashell <- runTopGoAnalysis(DEgeneSet = seashell$V1,
                                      dds = subA_dds2,
                                      GOdb = GOList,
                                      onts = "BP",
                                      nodeSize = 50)


top50_enriched_seashell <- exportGOtable(GO_seashell, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_seashell, "top50_enriched_seashell.tsv")

#pink
GO_pink <- runTopGoAnalysis(DEgeneSet = pink$V1,
                                dds = subA_dds2,
                                GOdb = GOList,
                                onts = "BP",
                                nodeSize = 50)


top50_enriched_pink <- exportGOtable(GO_pink, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_pink, "top50_enriched_pink.tsv")

#plum4
GO_plum4 <- runTopGoAnalysis(DEgeneSet = plum4$V1,
                            dds = subA_dds2,
                            GOdb = GOList,
                            onts = "BP",
                            nodeSize = 50)


top50_enriched_plum4 <- exportGOtable(GO_plum4, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_plum4, "top50_enriched_plum4.tsv")


#yellowgreen
GO_yellowgreen <- runTopGoAnalysis(DEgeneSet = yellowgreen$V1,
                             dds = subA_dds2,
                             GOdb = GOList,
                             onts = "BP",
                             nodeSize = 50)


top50_enriched_yellowgreen <- exportGOtable(GO_yellowgreen, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

write_tsv(top50_enriched_yellowgreen, "top50_enriched_yellowgreen.tsv")

#grey
GO_grey <- runTopGoAnalysis(DEgeneSet = grey$V1,
                                   dds = subA_dds2,
                                   GOdb = GOList,
                                   onts = "BP",
                                   nodeSize = 50)


top50_enriched_grey <- exportGOtable(GO_grey, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

head(top50_enriched_grey)

write_tsv(top50_enriched_grey, "top50_enriched_grey.tsv")


#royalblue4
GO_royalblue4 <- runTopGoAnalysis(DEgeneSet = royalblue4$V1,
                                  dds = subA_dds2,
                                  GOdb = GOList,
                                  onts = "BP",
                                  nodeSize = 50)


top50_enriched_royalblue4 <- exportGOtable(GO_royalblue4, "Fisher.weight01.unadjusted.p-value", n = 50) %>%
  mutate(Fisher.weight01 = as.numeric(Fisher.weight01))

head(top50_enriched_royalblue4)

write_tsv(top50_enriched_royalblue4, "top50_enriched_royalblue4.tsv")




###now let's make some figures so we don't spend a billion words on describing the GO terms. 

#subsets (top 20) GO terms for each comparison UP and DOWN, and lightpink3 genes

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/Pop4_bloomtime_PAPER/")

June12_top20 <- read.table("Jun12_subset_for_GO_graphing.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
June29_top20 <- read.table("Jun29_subset_for_GO_graphing.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
July18_top20 <- read.table("Jul18_subset_for_GO_graphing.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
hotpink_top20 <- read.table("hotpink_subset_for_GO_graphing.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
darkgreen_top20 <- read.table("darkgreen_subset_for_GO_graphing.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
gold_top20 <- read.table("gold_subset_for_GO_graphing.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
lavender_top20 <- read.table("lavender_subset_for_GO_graphing.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)

par(mfrow=c(2,1))

library(stringr)
June12_top20$Term2 <- str_wrap(June12_top20$Term, width = 40)

ggplot((June12_top20 %>%
         dplyr::filter(Direction=="UP")),aes(x=Term2,y=log.1.pvalue., fill=Direction)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6) +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,27),breaks=(c(0,5,10,15,20,25)))+
  scale_fill_manual(values=c("steelblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 5, 5))

ggplot((June12_top20 %>%
          dplyr::filter(Direction=="DOWN")),aes(x=Term2,y=log.1.pvalue., fill=Direction)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6) +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,27),breaks=(c(0,5,10,15,20,25)))+
  scale_fill_manual(values=c("pink1")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 5, 5))

###### June 29

June29_top20$Term2 <- str_wrap(June29_top20$Term, width = 40)

ggplot((June29_top20 %>%
          dplyr::filter(Direction=="UP")),aes(x=Term2,y=log.1.pvalue., fill=Direction)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6) +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,27),breaks=(c(0,5,10,15,20,25)))+
  scale_fill_manual(values=c("steelblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 5, 5))

ggplot((June29_top20 %>%
          dplyr::filter(Direction=="DOWN")),aes(x=Term2,y=log.1.pvalue., fill=Direction)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6) +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,27),breaks=(c(0,5,10,15,20,25)))+
  scale_fill_manual(values=c("pink1")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 5, 5))

###### July18

July18_top20$Term2 <- str_wrap(July18_top20$Term, width = 40)

ggplot((July18_top20 %>%
          dplyr::filter(Direction=="UP")),aes(x=Term2,y=log.1.pvalue., fill=Direction)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6) +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,32),breaks=(c(0,5,10,15,20,25)))+
  scale_fill_manual(values=c("steelblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 5, 5))

ggplot((July18_top20 %>%
          dplyr::filter(Direction=="DOWN")),aes(x=Term2,y=log.1.pvalue., fill=Direction)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6) +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,32),breaks=(c(0,5,10,15,20,25)))+
  scale_fill_manual(values=c("pink1")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 5, 5))

#### lightpink3 module GO terms
#install.packages("remotes")
library(ggpattern)

hotpink_top20$Term2 <- str_wrap(hotpink_top20$Term, width = 40)


ggplot(hotpink_top20,aes(x=Term2,y=log.1.pvalue.)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6, fill="hotpink", color="black") +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,27),breaks=(c(0,5,10,15,20,25)))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 5, 5))

darkgreen_top20$Term2 <- str_wrap(darkgreen_top20$Term, width = 40)


ggplot(darkgreen_top20,aes(x=Term2,y=log.1.pvalue.)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6, fill="darkgreen", color="black") +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,27),breaks=(c(0,5,10,15,20,25)))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 5, 5))

gold_top20$Term2 <- str_wrap(gold_top20$Term, width = 40)


ggplot(gold_top20,aes(x=Term2,y=log.1.pvalue.)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6, fill="gold", color="black") +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,27),breaks=(c(0,5,10,15,20,25)))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 5, 5))


lavender_top20$Term2 <- str_wrap(lavender_top20$Term, width = 37)

ggplot(lavender_top20,aes(x=Term2,y=log.1.pvalue.)) +
  geom_bar(stat="identity", show.legend = FALSE, width=0.6, fill="lavender", color="black") +
  geom_text(aes(label=Significant), vjust = -2.0, fontface=2, size=6) +
  geom_text(aes(label=round(Expected)), vjust = -0.5, size=5)+
  ylab("Log(1/p-value)") +
  xlab("") +
  scale_y_continuous(limits=c(0,27),breaks=(c(0,5,10,15,20,25)))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=18, vjust=0.5), axis.title.x=element_text(face="bold", size=0, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=20, margin=margin(r = 10)), axis.text.y = element_text(size=18), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"), plot.margin = margin(5, 5, 20, 5))




