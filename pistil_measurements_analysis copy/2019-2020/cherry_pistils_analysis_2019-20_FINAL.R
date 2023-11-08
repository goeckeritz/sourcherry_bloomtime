library(ggplot2)
library(tidyverse)
library(rmarkdown)
library(dplyr)

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/MicroscopeMeasurements/2019-2020/csvs/")

#Only have 4 columns here, so it's pretty basic. I'd like to add GDDs as a numerical and factor column, 
#bloom group (early or late), and we still need to normalize all trees to their individual pistil lengths. 
#note -- I must have forgotten to add 02-65 as 'done' for the last two dates, 5/11 and 5/15, because I collected open flowers on the 5/4. It's now been added

#MEAN Final pistil lengths:

early1 = read.csv("raw_FINAL_27-02-65_data.csv", header = T)
early1_final_length = mean(early1$Length..µm.)
early2 = read.csv("raw_FINAL_27-03-25_data.csv", header = T)
early2_final_length = mean(early2$Length..µm.)
early3 = read.csv("raw_FINAL_27-03-27_data.csv", header = T)
early3_final_length = mean(early3$Length)
early4 = read.csv("raw_FINAL_27-03-28_data.csv", header = T)
early4_final_length = mean(early4$Length..µm.)
early5 = read.csv("raw_FINAL_27-03-46_data.csv", header = T)
early5_final_length = mean(early5$Length..µm.)
early6 = read.csv("raw_FINAL_27-04-34_data.csv", header = T)
early6_final_length = mean(early6$Length..µm.)

late1 = read.csv("raw_FINAL_27-02-08_data.csv", header = T)
late1_final_length = mean(late1$Length..µm.)
late2 = read.csv("raw_FINAL_27-02-19_data.csv", header = T)
late2_final_length = mean(late2$Length..µm.)
late3 = read.csv("raw_FINAL_27-03-08_data.csv", header = T)
late3_final_length = mean(late3$Length..µm.)
late4 = read.csv("raw_FINAL_27-04-12_data.csv", header = T)
late4_final_length = mean(late4$Length..µm.)

#Combine all data, then create the pistil completion column in the grand file. 
temp = list.files(pattern = "raw_[0-9].*.csv")
temp
data = lapply(temp, read.csv)
str(data)
colnames = c("tree", "date", "Source", "Length", "Comment")
data = lapply(data, setNames, colnames)
library(dplyr)
pistils2020 = bind_rows(data)
library(tidyverse)
# we want columns "tree", "date", "Source", "Length", "Comment", "GDDs_since_Jan1_2019", "pistil_completion", "k", "j", 'bloom' to see if we can separate out its effects.)
x = vector(mode="numeric", length=4974)
additions = data.frame(matrix(ncol = 3, x))
additional_columns = c("GDDs_since_Jan1_2019", "pistil_completion", 'bloom')
colnames(additions) = additional_columns
pistils2020 = bind_cols(pistils2020, additions)

#divide each pistil measurement by the final length average for each tree! 
for (i in 1:nrow(pistils2020)){
  if (grepl("27-02-65", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/early1_final_length
  } else if (grepl("27-03-25", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/early2_final_length
  } else if (grepl("27-03-27", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/early3_final_length
  } else if (grepl("27-03-28", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/early4_final_length
  } else if (grepl("27-03-46", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/early5_final_length
  } else if (grepl("27-04-34", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/early6_final_length
  } else if (grepl("27-02-08", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/late1_final_length
  } else if (grepl("27-02-19", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/late2_final_length
  } else if (grepl("27-03-08", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/late3_final_length
  } else if (grepl("27-04-12", pistils2020$tree[i])){
    pistils2020[i, "pistil_completion"] <- pistils2020$Length[i]/late4_final_length
  }} 

#Fill in the GDDs for each collection date
for (i in 1:nrow(pistils2020)){
  if (grepl("10_13_2019", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2357
  } else if (grepl("11_16_19", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2423
  } else if (grepl("12_20_19", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2434
  } else if (grepl("2_22_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2456
  } else if (grepl("3_9_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2469
  } else if (grepl("3_20_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2483
  } else if (grepl("3_31_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2499
  } else if (grepl("4_6_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2516
  } else if (grepl("4_11_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2538
  } else if (grepl("4_18_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2551
  } else if (grepl("4_24_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2564
  } else if (grepl("4_27_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2579
  } else if (grepl("4_30_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2600
  } else if (grepl("5_2_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2616
  } else if (grepl("5_5_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2637
  } else if (grepl("5_11_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2653
  } else if (grepl("5_15_20", pistils2020$date[i])){
    pistils2020[i, "GDDs_since_Jan1_2019"] <- 2675
  }} 
#The GDDs should be classed as ordered factors - we do this because I want to know what time points early and late-blooming trees have statistical differences. 

pistils2020$GDDs_since_Jan1_2019_factor = ordered(as.factor(pistils2020$GDDs_since_Jan1_2019)) 
levels(pistils2020$GDDs_since_Jan1_2019_factor)
str(pistils2020)
  
#When the early-bloomers were finished blooming, there is an NA in the pistil_completion column which should be approximately 1
for (i in 1:nrow(pistils2020)){
  if (grepl("approximation", pistils2020$Source[i])){
    pistils2020[i, "pistil_completion"] <- 1
}}

for (i in 1:nrow(pistils2020)){
  if (grepl("27-02-65", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "early"
  } else if (grepl("27-03-25", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "early"
  } else if (grepl("27-03-27", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "early"
  } else if (grepl("27-03-28", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "early"
  } else if (grepl("27-03-46", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "early"
  } else if (grepl("27-04-34", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "early"
  } else if (grepl("27-02-08", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "late"
  } else if (grepl("27-02-19", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "late"
  } else if (grepl("27-03-08", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "late"
  } else if (grepl("27-04-12", pistils2020$tree[i])){
    pistils2020[i, "bloom"] <- "late"
  }} 

pistils2020$bloom = as.factor(pistils2020$bloom)
pistils2020$tree = as.factor(pistils2020$tree)

#Okay...I think we're ready to look at the data :D

library(fBasics)
library(lme4)
library(lmerTest)
library(multcomp)
library(carData)
library(car)
str(pistils2020)

#Our experimental design is the same as last year - but now we have more power with additional trees in each bloom group
#Let's first look at the residuals when we don't transform the data, shall we?

random_model1 <- lmer((pistil_completion) ~ bloom + GDDs_since_Jan1_2019_factor + bloom:GDDs_since_Jan1_2019_factor + (1|tree:bloom) + (1|tree:bloom:GDDs_since_Jan1_2019_factor), data=pistils2020)
fBasics::qqnormPlot(residuals(random_model1))
shapiro.test(residuals(random_model1))
hist(residuals(random_model1))
plot(residuals(random_model1)) #Lord... lol, fan shape indicates the need for a transformation. Let's do log transformation. 
#quite the narrow distribution with some outliers. Looks kinda normal actually. haha, but fan shape in residuals is BLEH

random_model2020 <- lmer((log(pistil_completion)) ~ bloom + GDDs_since_Jan1_2019_factor + bloom:GDDs_since_Jan1_2019_factor + (1|tree:bloom) + (1|tree:bloom:GDDs_since_Jan1_2019_factor), data=pistils2020)
fBasics::qqnormPlot(residuals(random_model2020)) #MUCH better
shapiro.test(residuals(random_model2020))
hist(residuals(random_model2020))
plot(residuals(random_model2020))
#Just like last year, the data is a bit skewed left but fairly normal looking. So we'll analyze it the similarly, but first drop extreme values (those below -0.5)
pistils2020$pistil_completion[which(residuals(random_model2020)<c(-0.5))] <- NA
#you'll have to rerun the script and recreate the pistils data frame fresh to get those values back. 
#Now we rerun the model. 
random_model2020 <- lmer((log(pistil_completion)) ~ bloom + GDDs_since_Jan1_2019_factor + bloom:GDDs_since_Jan1_2019_factor + (1|tree:bloom) + (1|tree:bloom:GDDs_since_Jan1_2019_factor), data=pistils2020)
fBasics::qqnormPlot(residuals(random_model2020)) #MUCH better
shapiro.test(residuals(random_model2020))
hist(residuals(random_model2020))
plot(residuals(random_model2020))
Anova(random_model2020, type="III")

#Time for mean separation of the fixed effects (particularly interested in the differences between bloom groups for each GDD level)
library(emmeans)

#This means we are comparing the means by collection/time point/GDD value. 
means2020 <- data.frame(emmeans(random_model2020, ~bloom + bloom:GDDs_since_Jan1_2019_factor))
str(means2020)

means2020 <- means2020[,c(1,2,3,6,7)] #drop SE and df columns for now
str(means2020)
means2020 = means2020 %>%
  mutate(GDDs_since_Jan1_2019_num = as.character(GDDs_since_Jan1_2019_factor)) %>%
  mutate(GDDs_since_Jan1_2019_num = as.numeric(GDDs_since_Jan1_2019_num))


#This is for the ggplot in a sec. #Though not really sure why this is needed. Lol
names(means2020)[1]<- c("Bloom")
means2020$Bloom <- as.character(means2020$Bloom)
means2020$Bloom[means2020$Bloom=="early"] <- "Early"
means2020$Bloom[means2020$Bloom=="late"] <- "Late"
str(means2020)


ggplot(means2020,aes(x=GDDs_since_Jan1_2019_factor,y=emmean)) +
  geom_point(aes(shape=Bloom, color=Bloom), size=5)+
  ggtitle("Carpel growth over time") + 
  xlab("2019 - 2020 Calendar Day") +
  ylab("Carpel Completion (log scale)") +
  scale_color_manual(values=c("black", "darkgrey"), labels = c("Early", "Late"), name = "Bloom Group") +
  scale_shape_manual(values=c(15,19), labels=c("Early", "Late"), name = "Bloom Group") +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL, color=Bloom),width=0.25) +
  scale_x_discrete(labels=c("2357" = "Oct 13", "2423" = "Nov 16", "2434" = "Dec 20", "2456" = "Feb 22",
                            "2469" = "Mar 9", "2483"= "Mar 20", "2499" = "Mar 31", 
                            "2516" = "Apr 6", "2538" = "Apr 11", "2551" = "Apr 18", "2564" = "Apr 24",
                            "2579" = "Apr 27", "2600" = "Apr 30", "2616" = "May 2", "2637" = "May 5", "2653" = "May 11", "2675" = "May 15")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))


# Problem with the above X scale is that it is just being plotted on a relative factor scale with no quantitative info

str(means2020)
str(pistils2020)
#stupid tricks to get the raw data on top of our averages:
levels(pistils2020$bloom)
pistils2020$log_Length = log(pistils2020$pistil_completion)
pistils2020$Bloom = pistils2020$bloom
pistils2020b = pistils2020 %>%
  dplyr::rename(emmean=log_Length) %>%
  dplyr::rename(GDDs_since_Jan1_2019_num = GDDs_since_Jan1_2019)
levels(pistils2020b$Bloom) = list("Early"="early", "Late"="late")
str(means2020)


ggplot(means2020,aes(x=GDDs_since_Jan1_2019_num,y=emmean)) +
  geom_point(aes(shape=Bloom, color=Bloom), size=4)+
  ggtitle("Carpel growth over time (2019-2020)") + 
  xlab("GDDs since January 1st 2019 (base temp 4.4C)") +
  ylab("Carpel Completion (log scale)") +
  xlim(2290,2800) +
  scale_color_manual(values=c("skyblue", "pink"), labels = c("Early", "Late"), name = "Bloom Group") +
  scale_shape_manual(values=c(15,19), labels=c("Early", "Late"), name = "Bloom Group") +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL, color=Bloom),width=0.25) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=15, margin=margin(r = 10)), axis.text.y = element_text(size=10), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15, margin=margin(b=10))) +
  geom_point(data=pistils2020b, aes(shape=Bloom, color=Bloom), size=0.75, position=position_jitter(width=2.5))


tests2020 <- data.frame(pairs(emmeans(random_model2020, ~bloom + bloom:GDDs_since_Jan1_2019_factor)))
pairs(emmeans(random_model2020, ~bloom + bloom:GDDs_since_Jan1_2019_factor)) #verified this is tukey adjusted
tests2020

str(tests2020)
head(tests2020)
tests2020 <- data.frame(tests2020, do.call(rbind,strsplit(as.character(tests2020$contrast), split="-")))
testsByGDD2020 <- list()
head(tests2020)
str(tests2020)

for(i in 1:nlevels(pistils2020$GDDs_since_Jan1_2019_factor)){
  testsByGDD2020[[i]] <- tests2020[grepl(
    as.character(levels(pistils2020$GDDs_since_Jan1_2019_factor))[i], tests2020$X1)&
      grepl(as.character(levels(pistils2020$GDDs_since_Jan1_2019_factor))[i], tests2020$X2), 1:6]
}

testsByGDD2020

#gotta make some dumb columns for plotting
str(pistils2020)
pistils2020$GDDs_Jan_num = pistils2020$GDDs_since_Jan1_2019
pistils2020$emmean = log(pistils2020$pistil_completion)
####Visualization by genotype - first, create a simple model. 
for_visualization <- lm(log(pistil_completion) ~ GDDs_since_Jan1_2019_factor + tree + tree:GDDs_since_Jan1_2019_factor, data=pistils2020)
means_tree <- data.frame(emmeans(for_visualization, ~tree + tree:GDDs_since_Jan1_2019_factor))
str(means_tree)
means_tree$GDDs_Jan_num <- as.numeric(as.character(means_tree$GDDs_since_Jan1_2019_factor))
str(means_tree)
levels(means_tree$tree)
means_tree$tree <- factor(means_tree$tree, levels = c('27-03-25', '27-03-46', '27-04-34','27-03-28', '27-02-65', '27-03-27', '27-03-08', '27-02-19', '27-04-12', '27-02-08'))
levels(means_tree$tree)

levels(pistils2020$tree)
#Rebecca thinks it would be nice to have the slopes on here. 
#I'll drop the pre-dormancy data points and just look at the growth from dormancy breakage to bloom. 
#That would be subsetting the data without 10/13 and 11/16.
#The last couple data points skew the line as well... ought to drop those too, I think. 

means_tree['bloom'] <- NA
for (i in 1:nrow(means_tree)){
  if (grepl("27-02-65", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "early"
  } else if (grepl("27-03-25", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "early"
  } else if (grepl("27-03-27", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "early"
  } else if (grepl("27-03-28", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "early"
  } else if (grepl("27-03-46", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "early"
  } else if (grepl("27-04-34", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "early"
  } else if (grepl("27-02-08", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "late"
  } else if (grepl("27-02-19", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "late"
  } else if (grepl("27-03-08", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "late"
  } else if (grepl("27-04-12", means_tree$tree[i])){
    means_tree[i, "bloom"] <- "late"
  }} 


#install.packages("ggnewscale")
library(ggnewscale)
  
  
ggplot(means_tree,aes(x=GDDs_Jan_num,y=emmean)) +
  geom_point(aes(shape=tree, col=tree), size=3)+
  ggtitle("Development Based on Pistil Growth (2019-2020)") + 
  xlab("GDDs since January 1st 2019 (base temp 4.4C)") +
  ylab("Pistil Completion (log scale)") +
  xlim(2290,2800) +
  scale_color_manual(values=c("steelblue", "steelblue2","steelblue1", "skyblue","lightblue2", "paleturquoise1", "pink2","violetred2", "violetred3", "violetred4"), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2", "27-04-34" = "Early 3", "27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-02-19" = "Late 2", "27-04-12" = "Late 3", "27-02-08"="Late 4"), name = "Bloom Group") +
  scale_shape_manual(values=c(15,15,15,15,15,15,19,19,19,19), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2", "27-04-34" = "Early 3", "27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-02-19" = "Late 2", "27-04-12" = "Late 3", "27-02-08"="Late 4"), name = "Bloom Group") +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL, color=tree),width=0.25) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15, face="bold"), 
        plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=15, face="bold"), legend.title=element_text(size=18, face="bold", margin=margin(b=10))) +
  geom_point(data=pistils2020, aes(shape=tree, color=tree), size=0.75, position=position_jitter(width=2.5)) +
  structure(ggplot2::standardise_aes_names("colour"), class = "new_aes") +
  # new_scale_color() +
  geom_smooth(data=subset(means_tree, GDDs_Jan_num > 2400 & GDDs_Jan_num < 2650), method='lm', se=FALSE, aes(group=bloom,color=bloom))+
  scale_color_manual(values = c("skyblue", "pink"), guide = NULL)

#we did it!!! yayyyyyyy 






