########################################
#                                      #
# Analyzing 2019-2020 winter ROS data  #
#                                      #
########################################

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/2019_winter_experiments/ROS_measurements/xylenol_orange/final_data/")
require(tidyverse)
require(ggplot2)

raw = read.csv("raw_data_BEST(60min_priority)min_10_22_edited.csv", header=TRUE)
str(raw)

#First, let's drop rows where we have -ive values (i.e., signal did not overcome background noise...or something odd happened that we can't pinpoint) and missing data (rows with NA)
raw = raw[raw$umoles_H2O2.gFW > 0, ] #12 wells with -ive values
raw = raw[!is.na(raw$umoles_H2O2.gFW),] #8 NAs - lost a rep for 04-12 and 03-28 

#Next, we need to get rid of outliers in the technical replicates within each tree within a plate. It might be cleanest to drop columns we aren't interested in for analysis first (like dilution factor and mg and such.)
names(raw)[1] = "genotype"
names(raw)[2] = "bloom_group"
names(raw)[4] = "number_of_buds"
names(raw)[18] = "umoles_H2O2_per_gFW"

raw_cleaned = subset(raw, select = -c(mg, time.min., equation, concentration_sample_well.umol.L., dilution_factor, concentration_sample.umol.L., sample_volume_extraction.L., notes))
str(raw_cleaned)
exists("collection", mode="variable")

library(fBasics)
library(lme4)
library(lmerTest)
library(broom)

#a measurement (technical replicate) is expected to be correlated with it's fellow technical reps - they are not representing any kind of biological variation.
#technically, am I dealing with a repeated measure here? If so, the design is similar to my forcing 
#experiments. Either way, a measurement will be nested within genotype at any given time point and the plate it's on. Meaning we'd have to create a dummy variable
#representing that combination. 

#seems I have two options here. Model these technical replicates as a random factor (lme4) or build a correlation matrix (nlme)
raw_cleaned$collection = as.factor(raw_cleaned$collection)
raw_cleaned$genotype = as.factor(raw_cleaned$genotype)
raw_cleaned$bloom_group = as.factor(raw_cleaned$bloom_group)
raw_cleaned$rep.plate = as.factor(raw_cleaned$rep.plate)
raw_cleaned$pseudoreplicate = as.factor(raw_cleaned$pseudoreplicate)
raw_cleaned = within(raw_cleaned, dummy <- as.factor(paste(genotype,chill_portions_from_Oct1,rep.plate, sep='_')))
raw_cleaned$chill_factor = ordered(as.factor(raw_cleaned$chill_portions_from_Oct1))
str(raw_cleaned)
levels(raw_cleaned$chill_factor)

library(nlme)

#treating technical replicates as repeated measures
model_rep_measures = lme(umoles_H2O2_per_gFW ~ bloom_group*chill_factor,
             random = ~1|dummy,
             correlation=corSymm(form=~1|dummy),
             raw_cleaned)
summary(model_rep_measures)
boxplot(model_rep_measures$residuals ~ raw_cleaned$bloom_group) #appears variance is approx. homogenous
hist(residuals(model_rep_measures)) # ~ very narrow distribution -- with some clear outliers; those are wacky technical replicates, no doubt. 
3*sd(residuals(model_rep_measures))
#remove wacky technical reps. 
raw_cleaned$umoles_H2O2_per_gFW[which(residuals(model_rep_measures)<c(-0.057))] <- NA  
raw_cleaned$umoles_H2O2_per_gFW[which(residuals(model_rep_measures)>c(0.057))] <- NA
raw_cleaned = raw_cleaned[!is.na(raw_cleaned$umoles_H2O2_per_gFW),]
#rerun the model
model_rep_measures = lme(umoles_H2O2_per_gFW ~ bloom_group*chill_factor,
                         random = ~1|dummy,
                         correlation=corSymm(form=~1|dummy),
                         raw_cleaned)
summary(model_rep_measures)
boxplot(model_rep_measures$residuals ~ raw_cleaned$bloom_group) #appears variance is approx. homogenous
library(car)
Anova(model_rep_measures, type="III") #degrees of freedom look good. 
#the only thing I don't love is that I can't have the random factor of genotype easily incorporated.


#reset the data to restore dropped values for the other model.
raw_cleaned = subset(raw, select = -c(mg, time.min., equation, concentration_sample_well.umol.L., dilution_factor, concentration_sample.umol.L., sample_volume_extraction.L., notes))
str(raw_cleaned)
exists("collection", mode="variable")
raw_cleaned$collection = as.factor(raw_cleaned$collection)
raw_cleaned$genotype = as.factor(raw_cleaned$genotype)
raw_cleaned$bloom_group = as.factor(raw_cleaned$bloom_group)
raw_cleaned$rep.plate = as.factor(raw_cleaned$rep.plate)
raw_cleaned$pseudoreplicate = as.factor(raw_cleaned$pseudoreplicate)
raw_cleaned = within(raw_cleaned, dummy <- as.factor(paste(genotype,chill_portions_from_Oct1,rep.plate, sep='_')))
raw_cleaned$chill_factor = ordered(as.factor(raw_cleaned$chill_portions_from_Oct1))
str(raw_cleaned)
levels(raw_cleaned$chill_factor)


#Model with just random effects:
model_random = lmer(umoles_H2O2_per_gFW ~ bloom_group*chill_factor + (1|genotype:bloom_group) + (1|dummy), data=raw_cleaned)
summary(model_random)
boxplot(residuals(model_random) ~ raw_cleaned$bloom_group) #appears variance is approx. homogenous
hist(residuals(model_random)) # ~ distribution looks essentially identical to the repeated measures model. 
3*sd(residuals(model_random)) #the standard deviation is a little lower - so it seems we're dealing with a similar model but it isn't the exact same. 
#remove wacky technical reps. 
raw_cleaned$umoles_H2O2_per_gFW[which(residuals(model_random)<c(-0.0572))] <- NA  
raw_cleaned$umoles_H2O2_per_gFW[which(residuals(model_random)>c(0.0572))] <- NA
raw_cleaned = raw_cleaned[!is.na(raw_cleaned$umoles_H2O2_per_gFW),]
#rerun our random model
model_random = lmer(umoles_H2O2_per_gFW ~ bloom_group*chill_factor + (1|genotype:bloom_group) + (1|dummy), data=raw_cleaned)
summary(model_random)
boxplot(residuals(model_random) ~ raw_cleaned$bloom_group) #appears variance is approx. homogenous
hist(residuals(model_random)) # ~ distribution looks essentially identical to the repeated measures model. 
Anova(model_random, type="III") #ahhh. hmm. bloom group is a bit less significant then our other model; but the differences in the p-values aren't too bad.  
#this might be because we've better accounted for the genotype-specific variation in this model.
#But how to know which model is 'better'?
#one thing is for certain: that high Pr(>Chisq) for bloom*chill factor indicates the genotypes are responding more or less the same across time points.
#which I find.. interesting.
#makes me wonder if chill sensing (if ros content is indeed a proxy for it) doesn't really depend on 
#development as they enter dormancy. More research needed!

#I definitely want to explicitly partition the variation within genotype (a random effect). So my gut says to go with this second model... 


#For the purposes of plotting some data, let's do some calculations. 
#We'll move forward by visualizing a few means; however, let's first find out what chill dates are statistically different from one another. 
library(emmeans)
library(pbkrtest)

#main effect of chill only
testsROS_chill <- data.frame(pairs(emmeans(model_random, ~chill_factor)))
str(testsROS_chill)
head(testsROS_chill)
#so, according to these results, the comparison p values are:
#6.6 v 26.6 portions = 0.000109 (the umoles of H2O2 per gram of fresh weight is significantly higher at 26.6 portions - which is the deep endodormancy most assuredly; 6.6 portions, there were still leaves and the flower buds were still actively developing).
#6.6 v 43.7 potions = 0.041 (the umoles of H2O2 per gram of fresh weight is significantly lower at 6.6 chill portions compared to 43.7 chill portions)
#no statistical difference for 6.6 and 74.8 chill portions
#ditto for 6.6 and 93.4.
#26.6 vs 74.8 portions = 0.0000084 - significantly higher H2O2 at 26.6 portions. 
#26.6 vs 93.4 = 0.0000018 - significantly higher H2O2 at 26.6 portions. 
#43.7 vs 74.8 chill portions = 0.011 - significantly higher H2O2 at 43.7 chill portions. 
#43.7 vs 93.4 = 0.0036 - significantly higher H2O2 at 43.7 chill portions. 
#no statistical difference between 74.8 and 93.4.

#I'd like plot mean values and color them by bloom group; for that, I need to make a model where
#genotype is a fixed effect. No need to reset the data here - only 8 observations were dropped anyway. 

model_temp = lmer(umoles_H2O2_per_gFW ~ genotype*chill_factor + (1|dummy), data=raw_cleaned)
summary(model_temp)
means_ROS_genotypes <- data.frame(emmeans(model_temp, ~chill_factor*genotype))
head(means_ROS_genotypes)
means_ROS_genotypes$chill_numeric = as.numeric(as.character(means_ROS_genotypes$chill_factor))
str(means_ROS_genotypes)
means_ROS_genotypes$genotype <- factor(means_ROS_genotypes$genotype, levels = c('27-03-25', '27-03-46', '27-04-34', '27-03-28', '27-02-65', '27-03-27', '27-03-08', '27-02-19', '27-04-12', '27-02-08'))
levels(means_ROS_genotypes$genotype)

ggplot(means_ROS_genotypes,aes(x=chill_numeric,y=emmean)) +
  geom_point(aes(color=genotype, shape=genotype), size=5)+
  geom_line(aes(color=genotype), size=1)+
  ggtitle("ROS Content During Dormancy") + 
  xlab("Chill portions from Oct 1st") +
  ylab("umoles of H2O2 per gFW") +
  scale_color_manual(values=c("steelblue", "steelblue2","steelblue1","skyblue","lightblue2", "paleturquoise1", "pink2","violetred2", "violetred3","violetred4"), labels = c("27-03-25"="Early 1", "27-03-46"="Early 2", "27-04-34" = "Early 3", "27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-02-19" = "Late 2", "27-04-12" = "Late 3", "27-02-08" = "Late 4"), name = "Bloom Group") +
  scale_shape_manual(values=c(15,15,15,15,15,15,19,19,19,19), labels = c("27-03-25"="Early 1", "27-03-46"="Early 2", "27-04-34" = "Early 3", "27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-02-19" = "Late 2", "27-04-12" = "Late 3", "27-02-08" = "Late 4"), name = "Bloom Group") +
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL, color=genotype),width=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15, face="bold"), plot.title=element_text(size=20, face = "bold", hjust=0.5), legend.text=element_text(size=15, face="bold"), legend.title=element_text(size=18, face="bold", margin=margin(b=10)))








