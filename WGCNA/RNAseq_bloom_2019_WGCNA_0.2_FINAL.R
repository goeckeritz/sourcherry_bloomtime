###got the FPKM matrix with some handy regular expressions. See file:~/Desktop/Desktop\ -\ Charity’s\ MacBook\ Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/get_fpkm_for_WGCNA.txt
###Let's do some correlation networks to find out where our DEGs of interest fall in the grand scheme of the biology... :D
###
#this first part isn't in the WGCNA tutorials and it seems to REALLY be fucking me over since a lot of my variances are low. I get rid of like 90% of my genes!!!
#after talking with Ben, I was reassured that we have some pretty good groups here, even with the 0.2 COV filter. 
#Since I'd like to keep as many genes as possible, I think I will stick to this cutoff. 
#The thing Ben seemed to confused about was why I was so concerned about connectivity metrics. 
#It's a fair question; as long as I have robust groupings that show up again and again, with biological relevance, 
#that's what matters most. Which I 1000% do. Plus, I could mess with these groupings all day, widdling the amount of genes down. Ben assured
#me it's nothing like p-hacking or something because the whole point of this is to find biological relevance. 
#So I significantly cleaned up this script according to his feedback. 
#Ben also recommended I try clustering with early- and late-bloomers separately, and not filter on covariance at all, and see what happens. 
#That will be in a separate script if I do that. 

#BiocManager::install("GO.db")
#BiocManager::install("preprocessCore")
#install.packages("doParallel")
library(doParallel)
library(DescTools)
library(WGCNA)
library(flashClust)
library(scales)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(BiocParallel)

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/WGCNA/post Manfeld:Goeckeritz meeting/")

#read in fpkm matrix:
fpkm = read.table("RNAseq_bloom_2019_FPKM.table", sep="\t", header=TRUE, stringsAsFactors = FALSE)
head(fpkm)
#read in functional info (we'll use it later)
functional = read.table("functional_annotation_QTL_round2_gene_names_B.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#need Gene_ID column in this file to drop the -R[A,B], to match Gene_ID in the fpkm file. 
functional$Gene_ID = gsub("-R[A,B]", "", functional$Gene_ID)

colnames(fpkm) <- sub("^X","",perl = T, colnames(fpkm))
colnames(fpkm) <- gsub(".","_",fixed = TRUE, colnames(fpkm))
colnames(fpkm) <- gsub("_FPKM","",fixed = TRUE, colnames(fpkm))
head(fpkm[,1:6])

# Processing Input Data: Filter Lowly Expressed Genes and Log Transform -----------------------------

#Make the row names the gene id
row.names(fpkm) <- fpkm$Gene_ID
nrow(fpkm)

#Remove genes where all samples have an FPKM value < 5; gotta get rid of this to get rid of low COVs 
fpkm.flt5 <- fpkm[,-1, drop=FALSE][apply(fpkm[,-1],1,function(x){!all(x < 5)}),] #had to add drop=FALSE because it was dropping my row names. Jesus WTH. Doesn't do that in Genevieve's script.
head(fpkm.flt5[,1:6])  
nrow(fpkm.flt5) #welp... bye bye to lots of my genes I suppose. We need to get rid of low ones because they are going to skew COVs

#boxplot of data
boxplot(fpkm.flt5) 
fpkm.flt5 <- log(fpkm.flt5+1) # +1 to get rid of negative values  
boxplot(fpkm.flt5) #oh hey that's much better. :D

# Processing Input Data: Filter by the Coefficient of Variance ------------

#Below is an example. We will test the function on the first row of our filtered data.
#We have to coerce the data to be numeric prior to pass it to the CoefVar function.
CoefVar(as.numeric(fpkm.flt5[1,])) #note: The higher the coefficient of variation, the greater the level of dispersion around the mean. It is generally expressed as a percentage.
#makes sense the variation would go down after log transformation.

#Let's try out a few different filters for quality control, and make a decision what to go with based on the number of genes left...
fpkm.flt5.cov0.2 <- fpkm.flt5[apply(fpkm.flt5, 1, function(x){!any(CoefVar(as.numeric(x)) < 0.2)}),] #how many genes have a coefficient of variation above 0.20?
head(fpkm.flt5.cov0.2)

#Lets look at the COV of the first row now. Remember, each row is a gene, and each column is a different tissue sample/replicate
CoefVar(as.numeric(fpkm.flt5.cov0.2[1,]))

#We're going to combine the three thresholds into a list to more easily work with the data ####
#### EDIT - I DID THAT IN A DIFFERENT SCRIPT AND DECIDED TO JUST GO WITH 0.2
#### however, I was too lazy to change the code to look like I was not longer working with a list of dataframes XD
fpkm.flt5.cov <- list("flt5.cov0.2"=fpkm.flt5.cov0.2)

#How many genes are remaining?
lapply(fpkm.flt5.cov, nrow)

# WGCNA: Initial Sample Clustering for QC ---------------------------------

#load the WGCNA and flastclust packages package if you haven't already; now we have pulled out genes with similar COVs, I think. 

#Transpose data matrix; WGCNA likes samples as rows, and columns as genes
fpkm.flt5.cov <- lapply(fpkm.flt5.cov,function(x){t(x)})

head(fpkm.flt5.cov$flt5.cov0.2[,1:5])

#Flagging genes with too many missing values
gsg <- lapply(fpkm.flt5.cov,goodSamplesGenes)

#Check is any genes were flagged that should be removed; Is it all okay? Should get TRUE if okay.
lapply(gsg, function(x){x$allOK})

#Generate sample tree; calculates distances of gene expression across samples, or something
sampleTree_bloom <- lapply(fpkm.flt5.cov,function(x){flashClust(dist(x), method = "average")})

str(sampleTree_bloom)

#Plot the tree; mar = margins. mfrow = #rows, #columns, cex = text size; might give an error if figure margins are too big
#par(mar = c(5,5,3,3),mfrow=c(1,3),cex=c(1.5)) #this makes thing uninterpretable *facepalm*

plot(sampleTree_bloom$flt5.cov0.2, sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 1.5,
     main = names(sampleTree_bloom)[1])

#note that I decided to leave 27-04-34_6_12_19_R2 out of this analysis since it was a HUGE outlier in the
#DESeq2 dataset. 
#What I tried a CoefVar of 0.02 up above and the clustering is always the same. 
#so, it seems 27-03-46_6_29_19_R3 is a bit of an outlier when looking at genes with variable expression across the samples. 

#getting rid of the outlier; I really didn't want to have to think to hard of about this, so a lot of this code is gonna be repetitive:
fpkm.flt5 <- fpkm[,-1, drop=FALSE][apply(fpkm[,-1],1,function(x){!all(x < 5)}),] 
fpkm.flt5 = dplyr::select(fpkm.flt5,-"27_03_46_6_29_19_R3") #if you look at our PCA for subA mapping, this library is the big triangle, on Jun29, and it does cluster closer with late-bloomers, so there's some sense in removing this noise. 
fpkm.flt5 <- log(fpkm.flt5+1)
fpkm.flt5.cov0.2 <- fpkm.flt5[apply(fpkm.flt5, 1, function(x){!any(CoefVar(as.numeric(x)) < 0.2)}),] #how many genes have a coefficient of variation above 0.20?
fpkm.flt5.cov <- list("flt5.cov0.2"=fpkm.flt5.cov0.2)
fpkm.flt5.cov <- lapply(fpkm.flt5.cov,function(x){t(x)})
sampleTree_bloom <- lapply(fpkm.flt5.cov,function(x){flashClust(dist(x), method = "average")})

plot(sampleTree_bloom$flt5.cov0.2, sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 1.5,
     main = names(sampleTree_bloom)[1])

#clustering looks better now -- replicates are grouping with one another. 
#0.2: same deal with 03-25. Ah, fascinating -- all the late-blooming data for 6/12 and 6/29 cluster in one big group. 7/18 is in one big group for the late-bloomers, but a bit scattered among the other time points for 04-34 and 03-46. hmmm.

#This explained why it's a good idea to filter on variance: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html 
#See data analysis question 2: Probesets or genes may be filtered by mean expression or variance (or their robust analogs such as median and median absolute deviation, MAD) since low-expressed or non-varying genes usually represent noise. 
#######################

# WGCNA: Pick the Soft Threshold Power ------------------------------------

#Testing different powers. Goal: have power less than 30, R squared greater than 0.8. This will optimize the scale free topology. Has everything to do with the decay function, which is the rate of decay of being able to tell your gene apart from another.. I think?? Connectivity decays as the power increases. 
powers = c(c(1:25), seq(from = 26, to=40, by=2)) #soft threshold is weighted
powers
#https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/

#Call the network topology analysis function
sft <- lapply(fpkm.flt5.cov,function(x){
  pickSoftThreshold(x,powerVector = powers, networkType = "signed", RsquaredCut = 0.85, verbose=5)
})
str(sft) 
#I guess the correlation matrix is raised to a power to reduce noise. inflate real signal, 'fuzz out'  
#the balance between reducing noise and... finding a good model fit? and ensuring it is a scale-free network? And this network's degrees should follow a power law.. the degree of a node in a network is the connections / edges it has to other nodes. sometimes those will be classified further into nodes going into the node of interest, and out of the node of interest. Most nodes are right-skewed, meaning that most nodes have a low # of connections, but a few have many - these are called hubs. 
#https://support.bioconductor.org/p/87024/#:~:text=The%20soft%20thresholding%2C%20is%20a,correlations%20in%20the%20adjacency%20matrix.:
#The soft thresholding, is a value used to power the correlation of the genes to that threshold. The assumption on that by raising the correlation to a power will reduce the noise of the correlations in the adjacency matrix.
#correlation -> adjacency (cor**power)-> Topological Overlap Measure -> clusters (with the DynamicTree algorithm). here is a hard threshold on the DynamicTree process, but it is the number of genes involved on each module.
#TOM calculation counts neighbors using a weighted sum: the weaker the connection, the less it counts.
#ah, well this is relatable: (For photography) A soft threshold is a preprocessing tool that reduces the BackGround in an image, so VoXels with intensity values below the threshold value are reduced (set to lower values, or even zero).

# Plot the results
par(mfrow = c(1,2),mar=c(5,5,3,3),cex=c(1.5))

# Scale-free topology fit index as a function of the soft-thresholding power
a <- plot(sft[[1]]$fitIndices[,1], -sign(sft[[1]]$fitIndices[,3])*sft[[1]]$fitIndices[,2],
          xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
          main = names(sft)[[1]],ylim = c(-0.5,1))

text(sft[[1]]$fitIndices[,1], -sign(sft[[1]]$fitIndices[,3])*sft[[1]]$fitIndices[,2],labels=powers,col="skyblue")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="firebrick")

# Mean connectivity as a function of the soft-thresholding power
a <- plot(sft[[1]]$fitIndices[,1], sft[[1]]$fitIndices[,5],
          xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",ylim = c(0,100),
          main = names(sft)[[1]])+text(sft[[1]]$fitIndices[,1], sft[[1]]$fitIndices[,5], labels=powers,col="red")

#'In co-expression networks, the connectivity measures how correlated a gene is with all other network genes.'
#https://edu.isb-sib.ch/pluginfile.php/158/course/section/65/_01_SIB2016_wgcna.pdf
#connectivity shouldn't be too high otherwise you've probably got something crazy driving your gene correlations.. but don't want it to be too low as we also want to actually observe clusters. 

#Check the values for each COV threshold
sft[[1]]$fitIndices[14:25,] #This is 0.1 COV -- average connectivity seems to be ~16, and powerEst is 18.
sft[[1]]$powerEstimate 

# WGCNA: Adjacency Matrix -------------------------------------------------

#Calculate the adjacency matrix
adjacency = adjacency(fpkm.flt5.cov$flt5.cov0.2, power = 14,type="signed") #type is the type of network you want. signed preserves signs of your PCC value. 
head(adjacency[,1:6])

#Export to cytoscape 
w <- exportNetworkToCytoscape(adjacency,edgeFile = "sour_cherry_bloom_cov0.2_flt5.txt", nodeFile = "sour_cherry_bloom_cov0.2_flt5.txt",threshold = 0.5)
#this created the files we will need to view our network in cytoscape; threshold defines whether an edge actually exists between two genes. below 0.5 usually results in toooooo many connections and you can't easily distinguish between modules. 
str(w)

# WGCNA: Dissimilarity Measure --------------------------------------------

#Create the topological overlap matrix similarity measure 
TOM <- TOMsimilarity(adjacency,TOMType = "signed")
#gene IDs are gone but everything is in the same order. 
head(TOM[,1:6])

#Turn the TOM into a dissimilarity measure
dissTOM <- 1-TOM #ha, fuck you TOM. 

#Notice how the diagnal is now 0 instead of 1
head(dissTOM[,1:6]) #so now lower numbers = more similar. 

# WGCNA: Hierarchical Clustering and Branch Cutting -------------------------------------------

#Build hierarchically clustered gene tree
geneTree <- hclust(as.dist(dissTOM),method = "average")
str(geneTree)
pheatmap(geneTree$merge)

#Cut the tree into modules, testing different cut heights. 
#This defines what a module is. well have to look at the tree to get a sense of how many clusters there ought to be
#A rough degree of control over what it means to be a subcluster is implemented by the parameter deepSplit.
#Ben recommended to go with the defaults - and the default is just over 0.99 - which, when I was running through this initially, was the height I decided to go with. 
dynamicMods.default <- cutreeDynamic(dendro = geneTree,distM = dissTOM, deepSplit = 2,minClusterSize = 30)

#Put them together in a list - here we're only doing default but I didn't want to adjst the code from when I was trouble shooting lol
dynamicMods <- list("default"=dynamicMods.default)

#Let's look at how many genes are in each module
lapply(dynamicMods,table)

#Part of WGCNA is that it uses colors to label modules and for graphing purposes. Here we are going to change the number labels to colors.
customColorOrder = c("yellowgreen", "pink", "lavender", "chocolate2","darkgreen","chartreuse4","gold","firebrick3","hotpink","lightblue2","dodgerblue3","royalblue4","khaki","black","plum4","palevioletred1","seashell", "lightcoral", "violetred3", "white")
dynamicColors <- lapply(dynamicMods, function(x){labels2colors(labels=x, colorSeq=customColorOrder)})
lapply(dynamicColors,table) 

#We will plot the gene tree with the module colors underneath to see how the modules compare to the tree.
plotDendroAndColors(geneTree,dynamicColors[[1]],dendroLabels = F,main=names(dynamicColors)[[1]],addGuide = T,hang=0.03,guideHang = 0.05)

#hmmm, there could be some merging here fore sure. 

# WGCNA: Module Merging ---------------------------------------------------
#Merging of modules whose expression profiles are very similar
#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent
#to merge such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire
#modules, we calculate their eigengenes (MEs) and cluster them on their consensus correlation, that is the minimum
#correlation across the two sets:

# Calculate eigengenes so that we can look at how they group, and decide if some modules should be merged. 
MEList <- lapply(dynamicColors,function(x){moduleEigengenes(fpkm.flt5.cov$flt5.cov0.2,colors = x)})

MEs <- lapply(MEList,function(x){x$eigengenes})

# Calculate dissimilarity of module eigengenes
MEDiss <- lapply(MEs,function(x){1-cor(x)})

# Cluster module eigengenes
METree <- lapply(MEDiss,function(x){hclust(as.dist(x),method = "average")})

#Plot the module trees
par(mar=c(2,2,2,2),mfrow=c(1,1))
a <- plot(METree[[1]],main = names(METree)[[1]],xlab = "")+abline(h=0.4,col="red")+abline(h=0.3,col="green")+abline(h=0.25,col="purple")+abline(h=0.2,col="blue")
#can we merge firebrick3 and royalblue4?

#Merge the modules according to the specified height
merge <- lapply(dynamicColors,function(x){mergeCloseModules(fpkm.flt5.cov$flt5.cov0.2,x,cutHeight = 0.20)}) #I want to keep lightpink3 intact.
str(merge)

mergedColors <- lapply(merge,function(x){x$colors})

str(mergedColors)

# Eigengenes of the new merged modules
mergedMEs <- lapply(merge,function(x){x$newMEs})

str(mergedMEs)

plotDendroAndColors(geneTree, cbind(dynamicColors[[1]], mergedColors[[1]]),
                    c("Initial Module Groups", "Merge cutHeight0.20"),
                    dendroLabels = FALSE, hang = 0.03, main=names(mergedColors)[1],
                    addGuide = TRUE, guideHang = 0.05)

#Let's look at how many genes are in each module now
lapply(mergedColors,table) 


# WGCNA: TOM and MDS Plots ------------------------------------------------

nGenes <- ncol(dissTOM)

#Set the number of genes to subset; hell let's just do em all. This is just for quality control anyway.  
nSelect <- 2071

#Set the random seed for consistancy
set.seed(10)

# sample indices
select <- sample(nGenes, size = nSelect)

#Use the randomly selected indices to select those genes from the dissTOM object
selectTOM = dissTOM[select, select]

#What are the dimensions of the subsetted dissTOM object?
dim(selectTOM)
head(selectTOM)

#Subset the module color assignments
selectColors <- lapply(mergedColors,function(x){x[select]})

#Calculate the multi-dimensional scaling
cmd1=cmdscale(as.dist(selectTOM),2)

#Plot the MDS plot; MDS = multi-dimensional scaling; sorta like PCA
#par(mfrow=c(2,2),mar=c(5,5,5,5))
plot(cmd1, col=as.character(selectColors[[1]]), main=paste("MDS Plot",names(selectColors)[[4]],sep = ": "),
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

# WGCNA: Extraction of the Genes and Expression Data for Each Module -------------------------

#Save filtered and transformed FPKM matrix; hooray! Need to format some stuff and extract certain things for later. (see below)
write.table(fpkm.flt5.cov$flt5.cov0.2,file = "gene_atlas_log_flt5_cov0.2_fpkm_matrix.txt",sep = "\t",quote = F)

# Select module probes
geneIDs <- colnames(fpkm.flt5.cov$flt5.cov0.2)
head(geneIDs)

#Convert the expression matrix to data frame; formatting! It was a number matrix before, we need it to be a data frame
cov.expr.df <- data.frame(fpkm.flt5.cov$flt5.cov0.2,stringsAsFactors = F)
head(cov.expr.df)
dim(cov.expr.df)
colnames(cov.expr.df) <- sub("^X","",perl = T, colnames(cov.expr.df))
colnames(cov.expr.df) <- gsub(".","-",fixed = TRUE, colnames(cov.expr.df))
colnames(cov.expr.df) <- gsub("_FPKM","",fixed = TRUE, colnames(cov.expr.df))
head(cov.expr.df[,1:6])


#Transpose so the samples are the columns again
cov.expr.df <- as.data.frame(t(cov.expr.df))

#Add a gene_id column from the row.names to use later with the left_join function
cov.expr.df$gene_id <- row.names(cov.expr.df)

head(cov.expr.df[,1:5])

#Indicate what the current path is..

path <- "/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/WGCNA/post Manfeld:Goeckeritz meeting/"


#Specify what the base name you want to be for the new directories
newDirectoryName <- "gene_atlas_flt5_cov0.2_modules_cut"
cutheight <- c("_default")

#Create an empty vector to hold our new directory names
newDir <- c()

#Loop through the three cut heights and paste the directory names together. Then, create the new directories.
for (i in 1:length(cutheight)){
  newDir[i] <- paste(path,newDirectoryName,cutheight[i],sep = "")
  dir.create(newDir[i])
}

#What directories are present now?
list.dirs(recursive = F)

#What are the names of the new directories?
newDir

#Loading dplyr
library(dplyr)

# Generating the function
print_module_colors <- function(x, y) {
  # First we are setting our working directory based on the name of the
  # element we are working with
  if (grepl("default", names(y)[x])) {
    setwd(newDir[1])
  } 
  
  # Next we extract the module color names
  colors <- unique(y[[x]])
  
  # Next we loop through each of the color names, performing the following
  # operations on each color
  for (i in 1:length(colors)) {
    
    # First we extract all the gene_ids which have been assigned to this module
    # color
    modGenes <- geneIDs[y[[x]] == colors[i]]
    
    # Then we coerce this vector into a data frame and add the expression
    # abundances using the left join function
    modExpr <- left_join(data.frame(gene_id = modGenes, stringsAsFactors = F), 
                         cov.expr.df, by = "gene_id")
    
    # Next we find the minimum and maximum expression value for each gene and
    # put those values into a vector
    x.min <- apply(modExpr[, -1], 1, function(x) {
      x[which.min(x)]
    })
    x.max <- apply(modExpr[, -1], 1, function(x) {
      x[which.max(x)]
    })
    
    # We will then use the min and max values to normalize the expression
    # abundances to values between 0 and 1. This normalization just helps with
    # visualizing the trends.
    modExpr.scale <- (modExpr[, -1] - (x.min))/(x.max - x.min)
    
    # We will assign the row names to be the gene_ids
    row.names(modExpr.scale) <- modExpr[, 1]
    
    # We will generate file names for these 6 data objects we just created
    # (the gene_ids by themselves, the expression abundances, and the
    # normalized/scaled expression abundances). Notice we use the 'colors'
    # object to have each name correspond to the module data.
    output_name_gene <- paste(colors[i], "_module_genes.txt", sep = "")
    output_name_expr <- paste(colors[i], "_module_log_flt5_cov0.2_fpkm_matrix.txt", 
                              sep = "")
    output_name_scale <- paste(colors[i], "_module_log_flt5_cov0.2_fpkm_scale_matrix.txt", 
                               sep = "")
    
    # Then we use these names we just generated to save the data to file
    write(modGenes, file = output_name_gene, ncolumns = 1, append = FALSE)
    write.table(modExpr, file = output_name_expr, quote = F, sep = "\t", 
                row.names = F)
    write.table(modExpr.scale, file = output_name_scale, sep = "\t", quote = F)
    
    # We will transpose the norm expression data
    modExpr.scale.t <- t(modExpr.scale)
    
    # We will now plot the norm expression data.  We will first create a new png
    # file with a name that includes the color
    pngName <- paste(colors[i], "_trend_plot.png", sep = "")
    png(pngName, type = "cairo")
    
    # This is just setting the plotting window to only incorporate 1 plot;
    # similar to mfrow
    split.screen(c(1, 1))
    
    # We will plot the first gene_id to set up the plot
    plot(modExpr.scale.t[, 1], ylim = c(0, 1), xlab = "Sample", ylab = "Z-score", 
         type = "l", lwd = 0.5, col = "grey")
    
    # Then we will loop through the rest of the columns (i.e. the gene_ids) to
    # plot the rest of the norm expression data for the rest of the gene_ids. We
    # are telling the graphics to NOT make a new plot for each gene and to NOT
    # make new axes for each gene.
    for (i in 1:length(modExpr.scale.t[1, ])) {
      screen(1, new = FALSE)
      plot(modExpr.scale.t[, i], ylim = c(0, 1), type = "l", lwd = 0.5, 
           col = "grey", xaxt = "n", yaxt = "n", ylab = "", xlab = "", 
           main = "", bty = "n")
    }
    
    # Finally, we close the png file.
    dev.off()
    
  }
}

# Execute function on each cut height object
f <- lapply(seq_along(mergedColors),print_module_colors,y=mergedColors)

#Let's look at what files are present for the black module.
list.files(pattern=glob2rx("black*"))

############################################## WGCNA: Extraction of the Connectivity Measures ----------------
# There are four connectivity measures that are captured for each gene: total connectivity (i.e. the sum of the entire row in the adjacency matrix), intramodular connectivity (i.e the sum of only the values for the genes in the same module), intermodular connectivity (i.e. the sum of only the values for the genes NOT in the same module), and the connectivity difference between the intra- and inter-modular connectivity.
# These values are an important metric of the quality of the module assignments and a good indication if a module should truly be considered a module. For example, an average intermodular connectivity (kOut) that is larger than the intramodular connectivity (kWithin), leading to a negative connectivity difference (kDiff), would suggest that genes in the module have more connections to genes in other modules than to genes within the module, defying the principle of a co-expression module. Such a module should be considered for removal or merging with other similar modules by going back to the module merging step and changing the height we used to merge the modules.
# Although, talking to Ben has changed my viewpoint on these and relaxed my apprehension of leaving negative modules alone.  

mod.edges <- lapply(mergedColors, function(x) {
  intramodularConnectivity(adjacency, x, scaleByMax = F)
})
str(mod.edges)

#writing the intramodular connectivities into a file for the different cut heights 
w <- lapply(seq_along(mod.edges), function(x, y) {
  write.table(y[[x]], file = paste("gene_atlas_flt5_cov0.2_intramodular_connectivity_", 
                                   names(y)[x], ".txt", sep = ""), sep = "\t", quote = F)
}, y = mod.edges)

#adding the module assignments to mod.edges
for (i in 1:length(mergedColors)) {
  mod.edges[[i]]$Module <- mergedColors[[i]]
}

str(mod.edges)

#group the genes together by module
q <- lapply(mod.edges, function(x) {
  x %>% group_by(Module) %>% summarise(avg = mean(kDiff), ratio = (mean(kDiff))/(mean(kTotal)))
})
q
##############################################

# WGCNA: Finding Module Biological Significance ----------------------------------

str(mergedMEs)

#The row names of each data frame corresponds to each of our samples. 
row.names(mergedMEs$default)

#For plotting purposes we need to assign sample info  to each of the samples and then add these  categories to the module eigengene data frame. We will make a data frame with our sample names. 
category <- data.frame("Sample"=row.names(fpkm.flt5.cov$flt5.cov0.2),stringsAsFactors = F)
head(category)

#Now we can add another column which has sample info
category$Category <- c(rep("Late2_6/12/19",3),rep("Late2_6/29/19",2),rep("Late2_7/18/19",3),
                       rep("Late1_6/12/19",3),rep("Late1_6/29/19",2),rep("Late1_7/18/19",3),
                       rep("Early1_6/12/19",3),rep("Early1_6/29/19",2),rep("Early1_7/18/19",3),
                       rep("Early2_6/12/19",3),rep("Early2_6/29/19",1),rep("Early2_7/18/19",3),
                       rep("Late3_6/12/19",3),rep("Late3_6/29/19",2),rep("Late3_7/18/19",3),
                       rep("Early3_6/12/19",2),rep("Early3_6/29/19",2),rep("Early3_7/18/19",3))

#Let's look at this data frame
category
row.names(category)

#The module eigengene table is in the format where each column corresponds to a module (i.e. it is in the short table format). We need to change this to the long format for plotting purposes where there is one column with the module name and another column with the eigengene value. We also need to make a new column with our sample ids to plot against them and we need to add our sample categories we just created. The below loop does all of these things for us for each of the cut heights.

library(tidyr)

for (i in 1:length(mergedMEs)){
  mergedMEs[[i]]$Sample <- as.numeric(as.character(row.names(mergedMEs[[i]]))) 
  mergedMEs[[i]]$Category <- category$Category
  mergedMEs[[i]] <- tidyr::gather(mergedMEs[[i]],key="module",value = "ME",-Sample,-Category)
}

mergedMEs[[1]]$Sample = as.numeric(rep(1:46),(length(unique(mergedMEs[[1]]$module))))
#You can see now that there are four columns: 1 for the sample id, 1 for the organ category, 1 for the module name, and 1 for the eigengene value.
str(mergedMEs)
columns <- as.data.frame(str_split_fixed(mergedMEs[[1]]$Category, '_', 2))
colnames(columns) = c("Genotype", "Date")
d <- cbind(mergedMEs[[1]], columns)
str(d)
head(d)
#reps_averaged_MEs = d %>%
#  group_by(Category) %>%
#  summarise(meanME=mean(ME), n=n())

#Now we can plot these values using ggplot2 and RColorBrewer. I will not go into details about how everything is working in ggplot since you will officially learn ggplot2 next week.
library(ggplot2)

#Change the plotting them to desired font sizes, types, legend positions, etc.
b <- theme(axis.text = element_text(size=12,colour="black"),      
           legend.text = element_text(size = 15),legend.title = element_text(size = 0),
           legend.title.align = (0.8),axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
           strip.text = element_text(size = 13,face = "bold"),legend.position = "right",
           element_line(linewidth=1))

#Let's look at how many genes are in each module for default cut height 
table(mergedColors$default)


#Plot the module eigengene value by sample
ggplot(mergedMEs[[1]],aes(Sample,ME,group=1,col=Category))+
  geom_hline(yintercept = 0,alpha=0.45)+
  geom_line(col="black")+geom_point()+facet_wrap(~module)+theme_bw()+
  scale_color_manual(values=c("lightskyblue","dodgerblue4","navyblue", "springgreen1", "springgreen4", "darkgreen", "thistle", "orchid4", "darkorchid4", "pink","hotpink", "violetred", "khaki", "goldenrod", "orange4", "indianred","red3","firebrick4"))+
  labs(y="Module Eigengene",title=names(mergedMEs)[[1]])+b

#Plot the module eigengene value by date
ggplot(d,aes(Date,ME,group=Genotype))+
  geom_point(aes(color=Genotype, shape=Genotype), size=2, stroke=1, position=position_jitter(width=0.05, height=0.0)) +
  geom_smooth(aes(color=Genotype), method=loess, size=0.75)+
  facet_wrap(~module)+theme_bw()+
  scale_color_manual(values=c("steelblue","steelblue2","steelblue1","pink","violetred2", "violetred3"))+
  scale_shape_manual(values=c(0,0,0,1,1,1))+
  labs(y="Module Eigengene",title=names(mergedMEs)[[1]]) +
  ggtitle("Eigengene Value by Date") +
  theme(axis.text.x=element_text(size=12, angle=45, face="bold", color="grey27", margin=margin(t=15)),
  strip.text = element_text(face = "bold", size=12), legend.title=element_text(size=16, face="bold"), legend.text=element_text(size=12,face="bold"),
  plot.title=element_text(size=18, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=16, margin=margin(r=5), face="bold"),
  axis.title.x=element_text(size=16, margin=margin(r=3), face="bold"))
  
###### OLD CODE BELOW
setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/WGCNA/post Manfeld:Goeckeritz meeting/gene_atlas_flt5_cov0.2_modules_cut_default/")
lightpink3_genes <- read.table("lightpink3_module_genes.txt", header=FALSE, stringsAsFactors = FALSE)
colnames(lightpink3_genes) = c("Gene_ID")
head(lightpink3_genes)
head(functional)

lightpink3_genes_functions <- left_join(x=lightpink3_genes, y=functional, by="Gene_ID") %>%
  dplyr::select(Gene_ID, "Putative_Function"="PANTHER...Pfam")

str(lightpink3_genes_functions)

write_tsv(lightpink3_genes_functions, "lightpink3_genes_and_functions.tsv")








#look at heatmap of gene expression for any module you want. 
setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/WGCNA/post Manfeld:Goeckeritz meeting/gene_atlas_flt5_cov0.2_modules_cut_default/")
#I can see the cluster of genes I want. They're right there at the top in the non-scaled heat map. 
lightpink3_def = read.table("lightpink3_module_log_flt5_cov0.2_fpkm_matrix.txt", stringsAsFactors = FALSE, header=TRUE)
head(lightpink3_def)
lightpink3_def = subset(lightpink3_def, select = -c(gene_id))
pheatmap(lightpink3_def)

darkgreen_def = read.table("darkgreen_module_log_flt5_cov0.2_fpkm_matrix.txt", stringsAsFactors = FALSE, header = TRUE)
head(darkgreen_def)
darkgreen_def = subset(darkgreen_def, select = -c(gene_id))
pheatmap(darkgreen_def)

lavender_def = read.table("lavender_module_log_flt5_cov0.2_fpkm_matrix.txt", stringsAsFactors = FALSE, header=TRUE)
head(lavender_def)
lavender_def = subset(lavender_def, select = -c(gene_id))
pheatmap(lavender_def)
#geez 03-25, calm the fuck down. 

#I would like to plot the eigengene values by date, like Ben did with time. 
#And if I could figure out a statistical test of significance, that'd also be good.


##### Visualization of module expression 
# we subsetted the fpkm matrix by module earlier; we can start with those. 

#################Probably a bit more trouble than it's worth; don't feel like it adds much my interpretation of the data when just looking at the ME - at least for lightpink3
setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/RNAseq_summer2019_bloom_analysis/WGCNA/post Manfeld:Goeckeritz meeting/gene_atlas_flt5_cov0.2_modules_cut_default/")
lightpink3_fpkm <- read.table("lightpink3_module_log_flt5_cov0.2_fpkm_scale_matrix.txt", header=TRUE, stringsAsFactors = FALSE)
head(lightpink3_fpkm)
colnames(lightpink3_fpkm) <- sub("^X","",perl = T, colnames(lightpink3_fpkm))
colnames(lightpink3_fpkm) <- gsub("27_02_19_","Late1-",fixed = TRUE, colnames(lightpink3_fpkm))
colnames(lightpink3_fpkm) <- gsub("27_03_08_","Late2-",fixed = TRUE, colnames(lightpink3_fpkm))
colnames(lightpink3_fpkm) <- gsub("27_04_12_","Late3-",fixed = TRUE, colnames(lightpink3_fpkm))
colnames(lightpink3_fpkm) <- gsub("27_03_25_","Early1-",fixed = TRUE, colnames(lightpink3_fpkm))
colnames(lightpink3_fpkm) <- gsub("27_03_46_","Early2-",fixed = TRUE, colnames(lightpink3_fpkm))
colnames(lightpink3_fpkm) <- gsub("27_04_34_","Early3-",fixed = TRUE, colnames(lightpink3_fpkm))
colnames(lightpink3_fpkm) <- gsub("_R.+","",fixed = FALSE, colnames(lightpink3_fpkm))
head(lightpink3_fpkm)
t_lightpink3_fpkm <- (t(lightpink3_fpkm))
head(t_lightpink3_fpkm)


column_2_add <- as.data.frame(as.factor(row.names(t_lightpink3_fpkm)))
names(column_2_add) = c('genotype_date')
head(column_2_add)
t_lightpink3_fpkm = cbind(t_lightpink3_fpkm, column_2_add)
row.names(t_lightpink3_fpkm) = NULL
head(t_lightpink3_fpkm)
t_lightpink3_fpkm = as.data.frame(t_lightpink3_fpkm)
gene_averages = aggregate(.~genotype_date, data=t_lightpink3_fpkm, mean)
more <- as.data.frame(str_split_fixed(gene_averages$genotype_date, '-', 2))
names(more) = c("genotype","date")
lightpink3_genes_to_plot <- cbind(gene_averages, more)
head(lightpink3_genes_to_plot)
lightpink3_gene_averages <- gather(lightpink3_genes_to_plot, key="Gene", value="avg_log2Fold_Change", 2:121)
head(lightpink3_gene_averages)
levels(lightpink3_gene_averages)

ggplot(lightpink3_gene_averages, aes(date, avg_log2Fold_Change))+
  geom_point(aes(color=genotype), size=0.5) +
  geom_line(aes(color=genotype, group=Gene))+
  scale_color_manual(values=c("steelblue","steelblue2","steelblue1","pink","violetred2", "violetred3"))+
  ggtitle("LightPink3 All Genes") +
  theme(axis.text.x=element_text(size=10, angle=45, face="bold", margin=margin(t=8)),
        strip.text = element_text(face = "bold", size = rel(1.0)), legend.title=element_text(size=14, face="bold"), 
        plot.title=element_text(size=14), axis.text.y=element_text(size=10), axis.title.y=element_text(size=14, margin=margin(r=5), face="bold"))







