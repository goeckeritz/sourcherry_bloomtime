#we back boiz, and it's time to look at the 2021-2022 pistil growth data. 

library(ggplot2)
library(tidyverse)
library(rmarkdown)
library(dplyr)
library(ggnewscale)

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/MicroscopeMeasurements/2021-2022/pistil_data_2021-2022/")

all_data = read.csv("combined_pistil_data_2021-2022b.csv", header=TRUE)
str(all_data)
levels(as.factor(all_data$Date))
final_lengths = read.csv("open_flower_averages.csv", header=TRUE)
str(final_lengths)

str(all_data)
all_data$GDDs_since_Jan1_2021_factor = ordered(as.factor(all_data$GDDs_since_Jan1_2021))
all_data$Genotype = as.factor(all_data$Genotype)
all_data$Length_uM = as.numeric(all_data$Length_uM) #the "done" was messing with the coding and making this column character. "done" should now be NAs
all_data$bloom = as.factor(all_data$bloom)

#time to normalize each value with final pistil length. 
for (i in 1:nrow(all_data)){
  if (grepl("27-02-08", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-02-08"])
  } else if (grepl("27-02-19", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-02-19"])
  } else if (grepl("27-02-51", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-02-51"])
  } else if (grepl("27-03-08", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-03-08"])  
  } else if  (grepl("27-03-24", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-03-24"])  
  } else if  (grepl("27-03-25", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-03-25"])  
  } else if  (grepl("27-03-27", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-03-27"])  
  } else if  (grepl("27-03-28", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-03-28"])  
  } else if  (grepl("27-03-46", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-03-46"])  
  } else if  (grepl("27-04-12", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-04-12"])  
  } else if  (grepl("27-04-40", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-04-40"])  
  } else if  (grepl("27-04-44", all_data$Genotype[i])){
    all_data[i, "normalized_Length"] <- all_data$Length_uM[i]/(final_lengths$Length_uM[final_lengths$Genotype=="27-04-44"])  
  }}
  
all_data$Genotype = gsub('-','_', all_data$Genotype)
head(all_data)
#when the tree has reached full bloom, it's 10 pistil values will be 1. 
all_data$normalized_Length = replace(all_data$normalized_Length, is.na(all_data$normalized_Length), 1)
head(all_data)  

#shall we take a look at the untransformed data?
library(emmeans)
library(fBasics)
library(lme4)
library(lmerTest)
library(car)
library(multcomp)

#simple model
model = lmer(normalized_Length ~ (1|Genotype:bloom:GDDs_since_Jan1_2021_factor) + (1|Genotype:bloom) + bloom + GDDs_since_Jan1_2021_factor + bloom:GDDs_since_Jan1_2021_factor, data=all_data)
plot(residuals(model))  #woof. fan-shaped as per usual.   

#simple log transformed model 
model_log = lmer(log(normalized_Length) ~ (1|Genotype:bloom:GDDs_since_Jan1_2021_factor) + (1|Genotype:bloom) + bloom + GDDs_since_Jan1_2021_factor + bloom:GDDs_since_Jan1_2021_factor, data=all_data)   
plot(residuals(model_log)) #much much better. 
hist(residuals(model_log)) #looks pretty damn normal, yay! But going to drop a few values as to not skew the model - I did this the last two years as well. 
all_data$normalized_Length[which(residuals(model_log)<c(-0.5))] <- NA  
all_data$normalized_Length[which(residuals(model_log)>c(0.5))] <- NA  


#rerun the model
model_log = lmer(log(normalized_Length) ~ (1|Genotype:bloom:GDDs_since_Jan1_2021_factor) + (1|Genotype:bloom) + bloom + GDDs_since_Jan1_2021_factor + bloom:GDDs_since_Jan1_2021_factor, data=all_data)   
plot(residuals(model_log))
hist(residuals(model_log))
fBasics::qqnormPlot(residuals(model_log))
Anova(model_log, type="III")

#we should only do anova with values of bloom within GDDs level. I still want to plot j means though. 
means = data.frame(emmeans(model_log, ~bloom + bloom*GDDs_since_Jan1_2021_factor))
means = means[,c(1,2,3,6,7)] #drop SE and df columns for now
#this next bit is adding the GDDs numeric column to our means data frame so we can plot the means against GDDs as a numeric value. 

means = means %>%
  mutate(GDDs_num = as.character(GDDs_since_Jan1_2021_factor)) %>%
  mutate(GDDs_num = as.numeric(GDDs_num))

str(means)
str(all_data)
levels(all_data$GDDs_since_Jan1_2021_factor)

levels(all_data$bloom)
all_data$log_Length = log(all_data$normalized_Length)
all_data$GDDs_num=as.numeric(all_data$GDDs_since_Jan1_2021)
##need x and y columns to match to lay the raw data over the means
all_data_b = all_data %>%
  dplyr::rename(emmean=log_Length) 
levels(means$bloom) = list("Early"="early", "Late"="late")
levels(all_data_b$bloom) = list("Early"="early", "Late"="late")
all_data_b$Genotype = as.factor(all_data_b$Genotype)
levels(all_data_b$Genotype)
all_data_b$Genotype = factor(gsub(pattern = "_", replacement = "-", x=all_data_b$Genotype))
str(all_data_b)
str(means)


ggplot(means,aes(x=GDDs_num,y=emmean)) +
  geom_point(aes(color=bloom, shape=bloom), size=4) +
  ggtitle("Carpel growth over time (2021-2022)") + 
  xlab("GDDs since January 1st 2021 (base temp 4.4C)") +
  ylab("Carpel Completion (log scale)") +
  scale_color_manual(values=c("skyblue", "pink"), labels = c("Early", "Late"), name = "Bloom Group") +
  scale_shape_manual(values=c(15,19), labels=c("Early", "Late"), name = "Bloom Group") +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL, color=bloom),width=0.25) +
  xlim(2850,3250) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))+
  geom_point(data=all_data_b, aes(shape=bloom, color=bloom), size=0.75, position=position_jitter(width=2.5))


#mean separation for bloom group. 
means_comp = data.frame(pairs(emmeans(model_log, ~bloom + bloom*GDDs_since_Jan1_2021_factor)))
means_comp <- data.frame(means_comp, do.call(rbind,strsplit(as.character(means_comp$contrast), split="-")))
means_comp_tests <- list()
head(means_comp)

for(i in 1:nlevels(all_data$GDDs_since_Jan1_2021_factor)){
  means_comp_tests[[i]] <- means_comp[grepl(
    as.character(levels(all_data$GDDs_since_Jan1_2021_factor))[i], means_comp$X1)&
      grepl(as.character(levels(all_data$GDDs_since_Jan1_2021_factor))[i], means_comp$X2), 1:6]
}

means_comp_tests 
write.csv(means_comp_tests, file="bloom_group_comparisons_by_GDDs_2021-22.csv")

#gotta make some dumb columns for plotting
str(all_data_b)
all_data_b$GDDs_Jan_num = all_data_b$GDDs_num
####Visualization by genotype - first, create a simple model. 
for_visualization <- lm(log(normalized_Length) ~ GDDs_since_Jan1_2021_factor + Genotype + Genotype:GDDs_since_Jan1_2021_factor, data=all_data_b)
means_Genotype <- data.frame(emmeans(for_visualization, ~Genotype + Genotype:GDDs_since_Jan1_2021_factor))
str(means_Genotype)
means_Genotype$GDDs_Jan_num <- as.numeric(as.character(means_Genotype$GDDs_since_Jan1_2021_factor))
str(means_Genotype)
levels(means_Genotype$Genotype)
means_Genotype$Genotype <- factor(means_Genotype$Genotype, levels = c('27-03-25', '27-03-46','27-03-28', '27-03-27', '27-03-24', '27-04-44', '27-02-51', '27-03-08', '27-02-19', '27-04-12', '27-02-08', '27-04-40'))
levels(means_Genotype$Genotype)
levels(all_data_b) <- factor(all_data_b$Genotype, levels = c('27-03-25', '27-03-46','27-03-28', '27-03-27', '27-03-24', '27-04-44', '27-02-51', '27-03-08', '27-02-19', '27-04-12', '27-02-08', '27-04-40'))
levels(all_data_b$Genotype)

means_Genotype['bloom'] <- NA
for (i in 1:nrow(means_Genotype)){
  if (grepl("27-03-25", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "early"
  } else if (grepl("27-03-27", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "early"
  } else if (grepl("27-03-28", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "early"
  } else if (grepl("27-03-46", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "early"
  } else if (grepl("27-02-08", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "late"
  } else if (grepl("27-02-19", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "late"
  } else if (grepl("27-03-08", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "late"
  } else if (grepl("27-04-12", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "late"
  } else if (grepl("27-04-44", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "early"
  } else if (grepl("27-02-51", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "early"
  } else if (grepl("27-04-40", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "late"
  } else if (grepl("27-03-24", means_Genotype$Genotype[i])){
    means_Genotype[i, "bloom"] <- "early"
  }} 


ggplot(means_Genotype,aes(x=GDDs_Jan_num,y=emmean)) +
  geom_point(aes(shape=Genotype, color=Genotype), size=3)+
  ggtitle("Development Based on Pistil Growth (2021-2022)") + 
  xlab("GDDs since January 1st 2021 (base temp 4.4C)") +
  ylab("Pistil Completion (log scale)") +
  #xlim(2750,3200) +
  scale_color_manual(values=c("steelblue", "steelblue2", "skyblue", "paleturquoise1", "lightsteelblue1", "lightblue3", "lightsteelblue4", "pink2","violetred2", "violetred3", "violetred4", "darkmagenta"), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2", "27-03-28" = "Early 4", "27-03-27" = "Early 6", "27-03-24"="Early8", "27-04-44" = "Early 9", "27-02-51"="Early 10", "27-03-08" = "Late 1", "27-02-19" = "Late 2", "27-04-12" = "Late 3", "27-02-08"="Late 4", "27-04-40"="Late 5"), name = "Bloom Group") +
  scale_shape_manual(values=c(15,15,15,15,15,15,15,19,19,19,19,19,19), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2", "27-03-28" = "Early 4", "27-03-27" = "Early 6", "27-03-24"="Early8", "27-04-44" = "Early 9", "27-02-51"="Early 10", "27-03-08" = "Late 1", "27-02-19" = "Late 2", "27-04-12" = "Late 3", "27-02-08"="Late 4", "27-04-40"="Late 5"), name = "Bloom Group") +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL, color=Genotype),width=0.25) +
  theme_minimal()+
  #geom_vline(xintercept=3184)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15, face="bold"), 
        plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=15, face="bold"), legend.title=element_text(size=18, face="bold", margin=margin(b=10))) +
  geom_point(data=all_data_b, aes(shape=Genotype, color=Genotype), size=0.75, position=position_jitter(width=2.5)) +
  structure(ggplot2::standardise_aes_names("colour"), class = "new_aes") +
  # new_scale_color() +
  geom_smooth(data=subset(means_Genotype, GDDs_Jan_num < 3185), method='lm', se=FALSE, aes(group=bloom,color=bloom))+
  scale_color_manual(values = c("skyblue", "pink"), guide = NULL)

#slope approximated for dormancy to the date when the first tree had bloomed.





