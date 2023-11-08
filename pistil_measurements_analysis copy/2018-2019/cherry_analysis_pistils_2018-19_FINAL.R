
library(ggplot2)
library(tidyverse)
library(rmarkdown)
library(dplyr)
library(ggnewscale)

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/MicroscopeMeasurements")
pistils <- read.csv("composite_data_2018-19_for_analysis.csv", header=TRUE)
head(pistils)

pistils <- pistils %>%
  dplyr::select(-Source) %>%
  dplyr::select(-Comment) %>%
  dplyr::rename(Date=Collection) %>%
  dplyr::rename(Days=Days_from_Jan_1st) #leave mids in; improves power

head(pistils)
str(pistils)
#copy GDDs from Jan 1st to a new column, name it GDDs_Jan_1st_numeric:
pistils$GDDs_Jan_1st_numeric = pistils$GDDs_from_Jan_1st_2018_C
#Then, replace it with a column where it is really GDDs from Sept20th - this is in case you want to plot the first time point as 0
pistils$GDDs_Sept20 = pistils$GDDs_Jan_1st_numeric - 2295
#We'll change the other original columns to factors for a factor levels effect model. 
pistils$GDDs_from_bloom_date_2018_C = ordered(as.factor(pistils$GDDs_from_bloom_date_2018_C))
levels(pistils$GDDs_from_bloom_date_2018_C)
pistils$GDDs_Jan = ordered(as.factor(pistils$GDDs_from_Jan_1st_2018_C)) #treating as factor so we can find out what dates developmental differences are significantly different. 
levels(pistils$GDDs_Jan)
pistils$Days = ordered(as.factor(pistils$Days))
str(pistils)
pistils$bloom = as.factor(pistils$bloom)
pistils$tree = as.factor(pistils$tree)
pistils$Date = as.factor(pistils$Date)



levels(pistils$Days) #Ah, there is one more level in Days because no GDDs happened between Feb21 and Mar 11. So GDDs of these two dates are equal and form one category in the GDDs_Jan model.  I'll check the GDDs again to see if there could even be 1 GDD difference; I would like to differentiate two different categories here. [checked, there isn't. So final analysis had the two dates meshed as one GDD category] 
levels(pistils$GDDs_Jan) #15 levels, Feb21 and Mar11 are the same level cuz they have the same amount of GDDs

library(fBasics)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(multcomp)
#pistil measurements are nested in genotype (tree) which is nested in early, late, and mid bloom categories (this year I had Balaton and Mont as part of the study - they are classified as 'mid')

random_model <- lmer(log(Length) ~ bloom + GDDs_Jan + bloom:GDDs_Jan + (1|tree:bloom) + (1|(tree:bloom):GDDs_Jan), data=pistils)
str(random_model)
fBasics::qqnormPlot(residuals(random_model))
shapiro.test(residuals(random_model))
hist(residuals(random_model))
#Drop extreme residuals
pistils$Length[which(residuals(random_model)<c(-0.5))] <- NA

#Re-run the model
random_model <- lmer(log(Length) ~ bloom + GDDs_Jan + bloom:GDDs_Jan + (1|tree:bloom) + (1|(tree:bloom):GDDs_Jan), data=pistils)
fBasics::qqnormPlot(residuals(random_model))
shapiro.test(residuals(random_model)) #Says it isn't normal, but we decided the residuals visually look normal, and that anovas are very robust and this should be a suitable analysis. Anova is reasonably robust against non-normality.
hist(residuals(random_model))

summary(random_model) #shows us the random effects
plot(residuals(random_model))
hist(residuals(random_model))

Anova(random_model, type="III") #Alllllll the significance, weeeeeee

#Now we do mean separation with the fixed effects. But we aren't worried about the genotype effects per se, more about the bloom time groups. 
#Note to self, add 1 to data where you didn't get measurements for mid and early groups (cuz they were finished blooming)
#This means we are comparing the means across all time points, or bloom categories within each GDDs_Jan level. 
meansBloom <- data.frame(emmeans(random_model, ~bloom + bloom:GDDs_Jan))
str(meansBloom)
#means <- emmeans(random_model, ~bloom + bloom:GDDs_Jan)
meansBloom <- meansBloom[,c(1,2,3,6,7)]
str(meansBloom)
names(meansBloom)[1]<- c("Bloom")
meansBloom$Bloom <- as.character(meansBloom$Bloom)
meansBloom$Bloom[meansBloom$Bloom=="early"] <- "Early"
meansBloom$Bloom[meansBloom$Bloom=="mid"] <- "Mid"
meansBloom$Bloom[meansBloom$Bloom=="late"] <- "Late"

no_mid <-meansBloom[meansBloom$Bloom!="Mid",]
ggplot(no_mid,aes(x=GDDs_Jan,y=emmean)) +
  geom_point(aes(shape=Bloom, color=Bloom), size=5)+
  ggtitle("Carpel growth over time") + 
  xlab("2018 - 2019 Calendar Day") +
  ylab("Carpel Completion (log scale)") +
  scale_color_manual(values=c("black", "darkgrey"), labels = c("Early", "Late"), name = "Bloom Group") +
  scale_shape_manual(values=c(15,19), labels=c("Early", "Late"), name = "Bloom Group") +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL, color=Bloom),width=0.25) +
  scale_x_discrete(labels=c("2295" = "Sep 20", "2538" = "Oct 30", "2558" = "Nov 28", "2559" = "Dec 27",
                            "2566" = "Feb 21 - Mar 11", "2572"= "Mar 18", "2588" = "Apr 1", 
                            "2628" = "Apr 15", "2672" = "Apr 23", "2691" = "Apr 26", "2711" = "May 4",
                            "2730" = "May 6", "2739" = "May 8", "2762" = "May 13", "2786" = "May 16")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

#Make comparable to the year 2020, but having the x axis give quantitative information on GDD accumulation. 
no_mid = no_mid %>%
  mutate(GDDs_Jan_num = as.character(GDDs_Jan)) %>%
  mutate(GDDs_Jan_num = as.numeric(GDDs_Jan_num)) %>%
  mutate(Bloom = as.factor(Bloom))

str(pistils)
levels(pistils$bloom)

#some dumb tricks to plot the raw data points on our averages. 
pistils_nomid = pistils %>%
  dplyr::filter(bloom!="mid")
levels(pistils_nomid$bloom)
pistils_nomid$log_Length = log(pistils_nomid$Length)
pistils_nomid = pistils_nomid %>%
  dplyr::rename(emmean=log_Length) %>%
  dplyr::rename(GDDs_Jan_num = GDDs_Jan_1st_numeric)

levels(pistils_nomid$bloom) = list("Early"="early", "Late"="late")
                        

str(no_mid)
str(pistils_nomid)

ggplot(no_mid,aes(x=GDDs_Jan_num,y=emmean)) +
  geom_point(aes(shape=Bloom, color=Bloom), size=4)+
  ggtitle("Carpel growth over time (2018-2019)") + 
  xlab("GDDs since January 1st 2018 (base temp 4.4C)") +
  ylab("Carpel Completion (log scale)") +
  xlim(2290,2800) +
  scale_color_manual(values=c("skyblue", "pink"), labels = c("Early", "Late"), name = "Bloom Group") +
  scale_shape_manual(values=c(15,19), labels=c("Early", "Late"), name = "Bloom Group") +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL, color=Bloom),width=0.25) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=15, margin=margin(r = 10)), axis.text.y = element_text(size=10), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15, margin=margin(b=10)))+
  geom_point(data=pistils_nomid, aes(shape=bloom, color=bloom), size=0.75, position=position_jitter(width=2.5))
  

#why you look funny, error bars??
  
#confidence_int <- confint(means)
#confidence_int
#Comparing the means within each time point, to find out when they are different. 

#How to do tests within each GDD, or Date or GDDs_Jan from Jan 1st to keep it ordered; pairs makes it a tukey test. 
tests <- data.frame(pairs(emmeans(random_model, ~bloom + bloom:GDDs_Jan)))
tests <- data.frame(tests, do.call(rbind,strsplit(as.character(tests$contrast), split="-")))
testsByGDD <- list()
head(tests)

for(i in 1:nlevels(pistils$GDDs_Jan)){
  testsByGDD[[i]] <- tests[grepl(
    as.character(levels(pistils$GDDs_Jan))[i], tests$X1)&
      grepl(as.character(levels(pistils$GDDs_Jan))[i], tests$X2), 1:6]
}

testsByGDD
  
####Visualization by genotype - first, create a simple model. 
for_visualization <- lm(log(Length) ~ GDDs_Jan + tree + tree:GDDs_Jan, data=pistils)
means_tree <- data.frame(emmeans(for_visualization, ~tree + tree:GDDs_Jan))
str(means_tree)
means_tree$GDDs_Jan_num <- as.numeric(as.character(means_tree$GDDs_Jan))
str(means_tree)
levels(means_tree$tree)
means_tree$tree <- factor(means_tree$tree, levels = c('27-03-25', '27-03-46', '27-04-34','27-03-08', '27-02-19', '27-04-12', 'mont','balaton'))
levels(means_tree$tree)


no_mid_trees <- means_tree[means_tree$tree!="mont",]
no_mid_trees <- no_mid_trees[no_mid_trees$tree!="balaton",]

no_mid_trees['bloom'] <- NA
for (i in 1:nrow(no_mid_trees)){
   if (grepl("27-03-25", no_mid_trees$tree[i])){
    no_mid_trees[i, "bloom"] <- "early"
  } else if (grepl("27-03-46", no_mid_trees$tree[i])){
    no_mid_trees[i, "bloom"] <- "early"
  } else if (grepl("27-04-34", no_mid_trees$tree[i])){
    no_mid_trees[i, "bloom"] <- "early"
    no_mid_trees[i, "bloom"] <- "late"
  } else if (grepl("27-02-19", no_mid_trees$tree[i])){
    no_mid_trees[i, "bloom"] <- "late"
  } else if (grepl("27-03-08", no_mid_trees$tree[i])){
    no_mid_trees[i, "bloom"] <- "late"
  } else if (grepl("27-04-12", no_mid_trees$tree[i])){
    no_mid_trees[i, "bloom"] <- "late"
  }} 



ggplot(no_mid_trees,aes(x=GDDs_Jan_num,y=emmean)) +
  geom_point(aes(shape=tree, color=tree), size=3)+
  ggtitle("Development Based on Pistil Growth (2018-2019)") + 
  xlab("GDDs since January 1st 2018 (base temp 4.4C)") +
  ylab("Pistil Completion (log scale)") +
  xlim(2290,2800) +
  scale_color_manual(values=c("steelblue", "steelblue2","steelblue1", "pink2","violetred2", "violetred3"), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2", "27-04-34" = "Early 3", "27-03-08" = "Late 1", "27-02-19" = "Late 2", "27-04-12" = "Late 3"), name = "Bloom Group") +
  scale_shape_manual(values=c(15,15,15,19,19,19), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2", "27-04-34" = "Early 3", "27-03-08" = "Late 1", "27-02-19" = "Late 2", "27-04-12" = "Late 3"), name = "Bloom Group") +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL, color=tree),width=0.25) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15, face="bold"), 
        plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=15, face="bold"), legend.title=element_text(size=18, face="bold", margin=margin(b=10))) +
  geom_point(data=pistils_nomid, aes(shape=tree, color=tree), size=0.75, position=position_jitter(width=2.5))+
  structure(ggplot2::standardise_aes_names("colour"), class = "new_aes") +
  # new_scale_color() +
  geom_smooth(data=subset(no_mid_trees, GDDs_Jan_num > 2558 & GDDs_Jan_num < 2750), method='lm', se=FALSE, aes(group=bloom,color=bloom))+
  scale_color_manual(values = c("skyblue", "pink"), guide = NULL)

#slope approximated for dormancy to the date when the first tree had bloomed.


