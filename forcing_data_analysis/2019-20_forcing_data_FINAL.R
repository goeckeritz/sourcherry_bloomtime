#we're starting fresh. We'll deal with repeated measures, but do several models where we look at stage 51+, 52+, 53+, etc, as a proportion.
#So, our response variable becomes continuous (thank God) Trying to do ordinal linear modeling and repeated measurements at the same time with clmm2 is basically impossible.

#This is the 2019-2020 forcing data. 

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/2019_winter_experiments/forcing_experiments/")

#using lme because it seems to be the only mixed linear effects modeling package that allows for a correlation matrix (repeated measures)
library(readxl)
library(ggplot2)
library(tidyverse)
library(nlme) #nlme allows you to specify var-cov matrices for random effects while lme4 does not. But people go back and forth on this in stack overflow O.o. We need to modify the VC structure for our repeated measures data!
library(emmeans)
library(phia)
library(car)
library(dplyr)
library(MASS)
#library(lme4)

#ultimately, we need to tell R that rep is nested within a genotype at each collection.
#to code that, I think we're gonna need to make a dummy variable that's a combination of genotype and chill_portions

data = read.csv("2019-2020_forcing_data_proportion.csv")
str(data)
#this data has been cleaned up in some ways (e.g., 04-34 was dropped cuz it behaved strangely and was hard to score due to bad flowers)

#Dropping observations/rows where I screwed up the count on accident 
data = data[(data$BBCH_50plus_proportion==1),]

str(data)
data = as.data.frame(data)

#I can drop a bunch of columns I don't need, then reclassify of some of the variables. 
#We'll model 'buds at stage 51 or above', 'buds at stage 52 or above', and so on, one model per stage. 
data = data %>%
  dplyr::select(-c('Hour', 'BBCH_50_proportion', 'BBCH_51_proportion', 'BBCH_52_proportion', 'BBCH_53_proportion', 'BBCH_54_proportion', 'BBCH_55_proportion', 'BBCH_56_proportion', 'BBCH_57_proportion', 'BBCH_59_proportion'))
data$chill_portions_since_oct1 = ordered(as.factor(data$chill_portions_since_oct1))
data$Genotype = as.factor(data$Genotype)
data$Group = as.factor(data$Group)
data$Date = as.factor(data$Date)
data$rep = ordered((as.factor(data$rep)),levels=c("R0","R2","R3", "R4", "R5","R6","R7","R8","R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16"))
data = within(data, dummy <- paste(Genotype,chill_portions_since_oct1, sep='_'))
data$dummy = as.factor(data$dummy) #we'll have to leave it unordered for simplicity. 
str(data)
levels(data$chill_portions_since_oct1) #order is correct.
is.ordered(data$chill_portions_since_oct1)
levels(data$rep)
isNested(data$rep, data$dummy) #we'll explicitly make it nested in our model.

#Now we need to figure out our model, and how to specify our var-cov matrix and which measures are repeated. 
#within each chill value, within each genotype/tree, scoring is repeated on every flower every couple of days.
#each repeated measurement is treated as a random effect and adjusted with the correlation matrix (to the best of my knowledge)

#for BBCH stage 51 (bud swell) and beyond! Note:Date is the repeated measure. 
model51= lme(BBCH_51plus_proportion ~ Group*chill_portions_since_oct1 + chill_portions_since_oct1*GDD_4.4C + Group*GDD_4.4C,
             random = ~1|dummy/rep,
             correlation=corSymm(form=~1|dummy/rep),
             data)

summary(model51)
boxplot(model51$residuals ~ data$Group) #appears variance is approx. homogenous
hist(residuals(model51)) # ~ normal, but a little left-skewed 
sd(residuals(model51))
Anova(model51, type='III')
#anova comparisons of models recommends having a chill*GDD and Group*chill AND Group*GDD interaction term.
#Group, chill, and GDDs are significantly affecting budbreak stage to 51 and beyond. Let's do separation of the marginal means of each group since interaction terms are not significant.
meansBBCH51 = data.frame(emmeans(model51, ~Group))
testsBBCH51 <- data.frame(pairs(emmeans(model51,  ~Group))) #104 is ~the average GDDs... would we call it the marginal mean of GDDs? working with numeric is a little hard to think about in this context. Mean separation seems a little weird here; let's just look at the means.
testsBBCH51 #no need to do my usual loop since we're only working with a single contrast. 
str(meansBBCH51)

ggplot(meansBBCH51,aes(x=Group,y=emmean, fill=Group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  ggtitle("Overall Mean Difference (2019-2020)") + 
  xlab("Bloom Group") +
  ylab("Stage 51+ Progression") +
  scale_fill_manual(values=c("steelblue", "lightpink")) +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=0.4)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=16), axis.title.x=element_text(face="bold", size=15, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=15, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.position="none")
#error bars are the upper and lower 95% confidence intervals for each bloom group. 
#if we really wanna, we can also plot the group means by chill collection:

meansBBCH51b = data.frame(emmeans(model51, ~Group*chill_portions_since_oct1))
str(meansBBCH51b)
meansBBCH51b = meansBBCH51b %>%
  mutate(chill_num = as.numeric(as.character(chill_portions_since_oct1)))
str(meansBBCH51b)

ggplot(meansBBCH51b,aes(x=chill_num,y=emmean)) +
  geom_point(aes(color=Group, shape=Group), size=6) +
  ggtitle("Average Stage 51+ Budbreak At Different Chill Levels (2019-2020)") +
  ylab("Progressed to Stage 51+ or further") +
  xlab("Chill Portions (Dynamic Model)") +
  scale_color_manual(values=c("steelblue","lightpink"), labels=c("Early","Late"), name="Bloom Group") +
  scale_shape_manual(values=c(15,19), labels=c("Early","Late"), name="Bloom Group")+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=2, color=Group)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=15, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5),
        legend.text=element_text(size=12), legend.title=element_text(size=15))

#it would be nice to be able to plot the raw data points for each bloom group, like I did with the pistils.
#so let's make a shitty model and do some tricks to do this. 

model51_genotypes= lme(BBCH_51plus_proportion ~ Genotype*chill_portions_since_oct1 + chill_portions_since_oct1*GDD_4.4C + Genotype*GDD_4.4C,
             random = ~1|dummy/rep,
             correlation=corSymm(form=~1|dummy/rep),
             data)
meansBBCH51_genotypes = data.frame(emmeans(model51_genotypes, ~Genotype*chill_portions_since_oct1)) #I suppose I don't need the confidence intervals here. 
str(meansBBCH51_genotypes)
meansBBCH51_genotypes = meansBBCH51_genotypes %>%
  mutate(chill_num = as.numeric(as.character(chill_portions_since_oct1)))
levels(meansBBCH51_genotypes$Genotype)
meansBBCH51_genotypes$Genotype <- factor(meansBBCH51_genotypes$Genotype, levels = c('27-03-25', '27-03-46', '27-03-28', '27-02-65', '27-03-27', '27-03-08', '27-04-12', '27-02-08'))
str(meansBBCH51_genotypes)
levels(meansBBCH51_genotypes$Genotype)
levels(meansBBCH51b$Group)
#wait a sec... I can simply combine these dataframes.
y = meansBBCH51b %>%
  dplyr::rename(Genotype=Group)
x <- rbind(meansBBCH51_genotypes,y)
levels(x$Genotype)

ggplot(x,aes(x=chill_num,y=emmean)) +
  geom_point(aes(color=Genotype, shape=Genotype, size=Genotype, stroke=0.75)) +
  labs(title="Average Stage 51+ Budbreak",
       subtitle="with Chill Accumulation (2019-2020)") +
  ylab("Progressed to Stage 51+ or further") +
  xlab("Chill Portions (Dynamic Model)") +
  scale_color_manual(values=c("steelblue", "steelblue2","skyblue","lightblue2", "paleturquoise1","pink2", "violetred3","violetred4", "steelblue","pink"), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group") +
  scale_shape_manual(values=c(12,12,12,12,12,10,10,10,15,19), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group")+  
  scale_size_manual(values=c(2,2,2,2,2,2,2,2,4,4), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group")+
  scale_x_continuous(breaks=c(20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=2, color=Genotype)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"))


################################## On to modeling stage 52

#In order to start fresh for the next model, the data needs to be read in again so the dropped observations return. 
data = read.csv("2019-2020_forcing_data_proportion.csv")
data = data[(data$BBCH_50plus_proportion==1),]
str(data)

data = data %>%
  dplyr::select(-c('Hour', 'BBCH_50_proportion', 'BBCH_51_proportion', 'BBCH_52_proportion', 'BBCH_53_proportion', 'BBCH_54_proportion', 'BBCH_55_proportion', 'BBCH_56_proportion', 'BBCH_57_proportion', 'BBCH_59_proportion'))
data$chill_portions_since_oct1 = ordered(as.factor(data$chill_portions_since_oct1))
data$Genotype = as.factor(data$Genotype)
data$Group = as.factor(data$Group)
data$Date = as.factor(data$Date)
data = within(data, dummy <- paste(Genotype,chill_portions_since_oct1, sep='_'))
data$dummy = as.factor(data$dummy) #we'll have to leave it unordered for simplicity. 
data$rep = ordered((as.factor(data$rep)),levels=c("R0","R2","R3", "R4", "R5","R6","R7","R8","R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16"))
str(data)
levels(data$chill_portions_since_oct1)

model52=lme(BBCH_52plus_proportion ~ Group*chill_portions_since_oct1 + chill_portions_since_oct1*GDD_4.4C + Group*GDD_4.4C,
             random = ~1|dummy/rep,
             correlation=corSymm(form=~1|dummy/rep),
             data)

summary(model52)
boxplot(model52$residuals ~ data$Group) #appears variance is approx. homogenous
hist(residuals(model52)) # ~ normal, but a little left-skewed 
plot(residuals(model52)) # that straight line of points essentially says the trees didn't make it to that stage.
Anova(model52, type="III")
#Group, GDDs, chill, and the interaction of chill and heat significantly affects a flower bud's ability to reach and pass stage 52. 
#Let's look at the overall mean difference between the early and late bloom groups.
#Then, we'll plot the means by chill, like we did for stage 51.

meansBBCH52 = data.frame(emmeans(model52, ~Group))
testsBBCH52 <- data.frame(pairs(emmeans(model52,  ~Group))) #104 is ~the average GDDs... would we call it the marginal mean of GDDs? working with numeric is a little hard to think about in this context. Mean separation seems a little weird here; let's just look at the means.
testsBBCH52 #no need to do my usual loop since we're only working with a single contrast. 
str(meansBBCH52)

ggplot(meansBBCH52,aes(x=Group,y=emmean, fill=Group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  ggtitle("Overall Mean Difference (2019-2020)") + 
  xlab("Bloom Group") +
  ylab("Stage 52+ Progression") +
  scale_fill_manual(values=c("steelblue", "lightpink")) +
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=0.4)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=16), axis.title.x=element_text(face="bold", size=15, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=15, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.position="none")
#error bars are the upper and lower 95% confidence intervals for each bloom group. 
#if we really wanna, we can also plot the group means by chill collection:

meansBBCH52b = data.frame(emmeans(model52, ~Group*chill_portions_since_oct1))
str(meansBBCH52b)
meansBBCH52b = meansBBCH52b %>%
  mutate(chill_num = as.numeric(as.character(chill_portions_since_oct1)))
str(meansBBCH52b)

ggplot(meansBBCH52b,aes(x=chill_num,y=emmean)) +
  geom_point(aes(color=Group, shape=Group), size=6) +
  ggtitle("Average Stage 52+ Budbreak At Different Chill Levels (2019-2020)") +
  ylab("Progressed to Stage 52+ or further") +
  xlab("Chill Portions (Dynamic Model)") +
  scale_color_manual(values=c("steelblue","lightpink"), labels=c("Early","Late"), name="Bloom Group") +
  scale_shape_manual(values=c(15,19), labels=c("Early","Late"), name="Bloom Group")+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=2, color=Group)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=15, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5),
        legend.text=element_text(size=12), legend.title=element_text(size=15))

#it would be nice to be able to plot the raw data points for each bloom group, like I did with the pistils.
#so let's make a shitty model and do some tricks to do this. 

model52_genotypes= lme(BBCH_52plus_proportion ~ Genotype*chill_portions_since_oct1 + chill_portions_since_oct1*GDD_4.4C + Genotype*GDD_4.4C,
                       random = ~1|dummy/rep,
                       correlation=corSymm(form=~1|dummy/rep),
                       data)
meansBBCH52_genotypes = data.frame(emmeans(model52_genotypes, ~Genotype*chill_portions_since_oct1)) #I suppose I don't need the confidence intervals here. 
str(meansBBCH52_genotypes)
meansBBCH52_genotypes = meansBBCH52_genotypes %>%
  mutate(chill_num = as.numeric(as.character(chill_portions_since_oct1)))
levels(meansBBCH52_genotypes)
meansBBCH52_genotypes$Genotype <- factor(meansBBCH52_genotypes$Genotype, levels = c('27-03-25', '27-03-46', '27-03-28', '27-02-65', '27-03-27', '27-03-08', '27-04-12', '27-02-08'))
str(meansBBCH52_genotypes)
levels(meansBBCH52_genotypes$Genotype)
levels(meansBBCH52b$Group)
#wait a sec... I can simply combine these dataframes.
y = meansBBCH52b %>%
  dplyr::rename(Genotype=Group)
x <- rbind(meansBBCH52_genotypes,y)
levels(x$Genotype)

ggplot(x,aes(x=chill_num,y=emmean)) +
  geom_point(aes(color=Genotype, shape=Genotype, size=Genotype, stroke=0.75)) +
  labs(title="Average Stage 52+ Budbreak",
       subtitle="with Chill Accumulation (2019-2020)") +
  ylab("Progressed to Stage 52+ or further") +
  xlab("Chill Portions (Dynamic Model)") +
  scale_color_manual(values=c("steelblue", "steelblue2","skyblue","lightblue2", "paleturquoise1","pink2", "violetred3","violetred4", "steelblue","pink"), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group") +
  scale_shape_manual(values=c(12,12,12,12,12,10,10,10,15,19), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group")+  
  scale_size_manual(values=c(2,2,2,2,2,2,2,2,4,4), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group")+
  scale_x_continuous(breaks=c(20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=2, color=Genotype)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"))

############# On to stage 53
data = read.csv("2019-2020_forcing_data_proportion.csv")
data = data[(data$BBCH_50plus_proportion==1),]
str(data)

data = data %>%
  dplyr::select(-c('Hour', 'BBCH_50_proportion', 'BBCH_51_proportion', 'BBCH_52_proportion', 'BBCH_53_proportion', 'BBCH_54_proportion', 'BBCH_55_proportion', 'BBCH_56_proportion', 'BBCH_57_proportion', 'BBCH_59_proportion'))
data$chill_portions_since_oct1 = ordered(as.factor(data$chill_portions_since_oct1))
data$Genotype = as.factor(data$Genotype)
data$Group = as.factor(data$Group)
data$Date = as.factor(data$Date)
data = within(data, dummy <- paste(Genotype,chill_portions_since_oct1, sep='_'))
data$dummy = as.factor(data$dummy) #we'll have to leave it unordered for simplicity. 
data$rep = ordered((as.factor(data$rep)),levels=c("R0","R2","R3", "R4", "R5","R6","R7","R8","R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16"))
str(data)

model53= lme(BBCH_53plus_proportion ~ Group*chill_portions_since_oct1 + chill_portions_since_oct1*GDD_4.4C + Group*GDD_4.4C,
             random = ~1|dummy/rep,
             correlation=corSymm(form=~1|dummy/rep),
             data)
summary(model53)
hist(residuals(model53))
plot(residuals(model53))
Anova(model53, type="III")

meansBBCH53 = data.frame(emmeans(model53, ~Group*chill_portions_since_oct1 + Group*GDD_4.4C))
str(meansBBCH53)
meansBBCH53 = meansBBCH53 %>%
  mutate(chill_num = as.numeric(as.character(chill_portions_since_oct1)))
str(meansBBCH53)

ggplot(meansBBCH53,aes(x=chill_num,y=emmean)) +
  geom_point(aes(color=Group, shape=Group), size=6) +
  ggtitle("Average Stage 53+ Budbreak At Different Chill Levels (2019-2020)") +
  ylab("Progressed to Stage 53+ or further") +
  xlab("Chill Portions (Dynamic Model)") +
  scale_color_manual(values=c("steelblue","lightpink"), labels=c("Early","Late"), name="Bloom Group") +
  scale_shape_manual(values=c(15,19), labels=c("Early","Late"), name="Bloom Group")+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=2, color=Group)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=15, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5),
        legend.text=element_text(size=12), legend.title=element_text(size=15))

testsBBCH53 <- data.frame(pairs(emmeans(model53, ~Group*chill_portions_since_oct1 + Group*GDD_4.4C)))

testsBBCH53 <- data.frame(testsBBCH53, do.call(rbind,strsplit(as.character(testsBBCH53$contrast), split="-")))
head(testsBBCH53)
str(testsBBCH53)
testsByChill53 <- list()

for(i in 1:nlevels(data$chill_portions_since_oct1)){
  testsByChill53[[i]] <- testsBBCH53[grepl(
    as.character(levels(data$chill_portions_since_oct1))[i], testsBBCH53$X1)&
      grepl(as.character(levels(data$chill_portions_since_oct1))[i], testsBBCH53$X2), 1:6]
}

testsByChill53 

##################plot with genotype averages on them:
model53_genotypes= lme(BBCH_53plus_proportion ~ Genotype*chill_portions_since_oct1 + chill_portions_since_oct1*GDD_4.4C + Genotype*GDD_4.4C,
                       random = ~1|dummy/rep,
                       correlation=corSymm(form=~1|dummy/rep),
                       data)
meansBBCH53_genotypes = data.frame(emmeans(model53_genotypes, ~Genotype*chill_portions_since_oct1)) #I suppose I don't need the confidence intervals here. 
str(meansBBCH53_genotypes)
meansBBCH53_genotypes = meansBBCH53_genotypes %>%
  mutate(chill_num = as.numeric(as.character(chill_portions_since_oct1)))
levels(meansBBCH53_genotypes)
meansBBCH53_genotypes$Genotype <- factor(meansBBCH53_genotypes$Genotype, levels = c('27-03-25', '27-03-46', '27-03-28', '27-02-65', '27-03-27', '27-03-08', '27-04-12', '27-02-08'))
str(meansBBCH53_genotypes)
levels(meansBBCH53_genotypes$Genotype)
levels(meansBBCH53b$Group)
#wait a sec... I can simply combine these dataframes.
y = meansBBCH53b %>%
  dplyr::rename(Genotype=Group)
x <- rbind(meansBBCH53_genotypes,y)
levels(x$Genotype)

ggplot(x,aes(x=chill_num,y=emmean)) +
  geom_point(aes(color=Genotype, shape=Genotype, size=Genotype, stroke=0.75)) +
  labs(title="Average Stage 53+ Budbreak",
       subtitle="with Chill Accumulation (2019-2020)") +
  ylab("Progressed to Stage 53+ or further") +
  xlab("Chill Portions (Dynamic Model)") +
  scale_color_manual(values=c("steelblue", "steelblue2","skyblue","lightblue2", "paleturquoise1","pink2", "violetred3","violetred4", "steelblue","pink"), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group") +
  scale_shape_manual(values=c(12,12,12,12,12,10,10,10,15,19), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group")+  
  scale_size_manual(values=c(2,2,2,2,2,2,2,2,4,4), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group")+
  scale_x_continuous(breaks=c(20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=2, color=Genotype)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"))


################################ on to stage 54!

data = read.csv("2019-2020_forcing_data_proportion.csv")
data = data[(data$BBCH_50plus_proportion==1),]
str(data)

data = data %>%
  dplyr::select(-c('Hour', 'BBCH_50_proportion', 'BBCH_51_proportion', 'BBCH_52_proportion', 'BBCH_53_proportion', 'BBCH_54_proportion', 'BBCH_55_proportion', 'BBCH_56_proportion', 'BBCH_57_proportion', 'BBCH_59_proportion'))
data$chill_portions_since_oct1 = ordered(as.factor(data$chill_portions_since_oct1))
data$Genotype = as.factor(data$Genotype)
data$Group = as.factor(data$Group)
data$Date = as.factor(data$Date)
data = within(data, dummy <- paste(Genotype,chill_portions_since_oct1, sep='_'))
data$dummy = as.factor(data$dummy) #we'll have to leave it unordered for simplicity. 
data$rep = ordered((as.factor(data$rep)),levels=c("R0","R2","R3", "R4", "R5","R6","R7","R8","R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16"))
str(data)

model54= lme(BBCH_54plus_proportion ~ Group*chill_portions_since_oct1 + chill_portions_since_oct1*GDD_4.4C + Group*GDD_4.4C,
             random = ~1|dummy/rep,
             correlation=corSymm(form=~1|dummy/rep),
             data)
summary(model54)
hist(residuals(model54))
plot(residuals(model54))
Anova(model54, type="III")

meansBBCH54 = data.frame(emmeans(model54, ~Group*chill_portions_since_oct1 + Group*GDD_4.4C)) #note, it doesn't really matter if we put chill*GDD here since we're looking at group means. 
str(meansBBCH54)
meansBBCH54 = meansBBCH54 %>%
  mutate(chill_num = as.numeric(as.character(chill_portions_since_oct1)))
str(meansBBCH54)

ggplot(meansBBCH54,aes(x=chill_num,y=emmean)) +
  geom_point(aes(color=Group, shape=Group), size=6) +
  ggtitle("Average Stage 54+ Budbreak At Different Chill Levels (2019-2020)") +
  ylab("Progressed to Stage 54+ or further") +
  xlab("Chill Portions (Dynamic Model)") +
  scale_color_manual(values=c("steelblue","lightpink"), labels=c("Early","Late"), name="Bloom Group") +
  scale_shape_manual(values=c(15,19), labels=c("Early","Late"), name="Bloom Group")+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=2, color=Group)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=15, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5),
        legend.text=element_text(size=12), legend.title=element_text(size=15))

testsBBCH54 <- data.frame(pairs(emmeans(model54, ~Group*chill_portions_since_oct1 + Group*GDD_4.4C)))

testsBBCH54 <- data.frame(testsBBCH54, do.call(rbind,strsplit(as.character(testsBBCH54$contrast), split="-")))
head(testsBBCH54)
str(testsBBCH54)
testsByChill54 <- list()

for(i in 1:nlevels(data$chill_portions_since_oct1)){
  testsByChill54[[i]] <- testsBBCH54[grepl(
    as.character(levels(data$chill_portions_since_oct1))[i], testsBBCH54$X1)&
      grepl(as.character(levels(data$chill_portions_since_oct1))[i], testsBBCH54$X2), 1:6]
}

testsByChill54

####to get plots with genotype averages:
model54_genotypes= lme(BBCH_54plus_proportion ~ Genotype*chill_portions_since_oct1 + chill_portions_since_oct1*GDD_4.4C + Genotype*GDD_4.4C,
                       random = ~1|dummy/rep,
                       correlation=corSymm(form=~1|dummy/rep),
                       data)
meansBBCH54_genotypes = data.frame(emmeans(model54_genotypes, ~Genotype*chill_portions_since_oct1*GDD_4.4C)) #I suppose I don't need the confidence intervals here. 
str(meansBBCH54_genotypes)
str(meansBBCH54)
meansBBCH54_genotypes = meansBBCH54_genotypes %>%
  mutate(chill_num = as.numeric(as.character(chill_portions_since_oct1)))
levels(meansBBCH54_genotypes)
meansBBCH54_genotypes$Genotype <- factor(meansBBCH54_genotypes$Genotype, levels = c('27-03-25', '27-03-46', '27-03-28', '27-02-65', '27-03-27', '27-03-08', '27-04-12', '27-02-08'))
str(meansBBCH54_genotypes)
levels(meansBBCH54_genotypes$Genotype)
levels(meansBBCH54$Group)
#wait a sec... I can simply combine these dataframes.
y = meansBBCH54 %>%
  dplyr::rename(Genotype=Group)
x <- rbind(meansBBCH54_genotypes,y)
levels(x$Genotype)

ggplot(x,aes(x=chill_num,y=emmean)) +
  geom_point(aes(color=Genotype, shape=Genotype, size=Genotype, stroke=0.75)) +
  labs(title="Average Stage 54+ Budbreak",
       subtitle="with Chill Accumulation (2019-2020)") +
  ylab("Progressed to Stage 54+ or further") +
  xlab("Chill Portions (Dynamic Model)") +
  scale_color_manual(values=c("steelblue", "steelblue2","skyblue","lightblue2", "paleturquoise1","pink2", "violetred3","violetred4", "steelblue","pink"), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group") +
  scale_shape_manual(values=c(12,12,12,12,12,10,10,10,15,19), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group")+  
  scale_size_manual(values=c(2,2,2,2,2,2,2,2,4,4), labels=c("27-03-25"="Early 1", "27-03-46"="Early 2","27-03-28" = "Early 4", "27-02-65" = "Early 5", "27-03-27" = "Early 6", "27-03-08" = "Late 1", "27-04-12" = "Late 3", "27-02-08" = "Late 4", "early"="Early Group Avg", "late"="Late Group Avg"), name="Bloom Group")+
  scale_x_continuous(breaks=c(20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL,width=2, color=Genotype)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15), axis.title.x=element_text(face="bold", size=18, margin=margin(t = 10)), axis.title.y=element_text(face="bold", size=18, margin=margin(r = 10)), axis.text.y = element_text(size=15), plot.title=element_text(size=22, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=18, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.title=element_text(size=18, face="bold"))




################## the rest of the script is simply examination of the raw data.

data = read.csv("2019-2020_forcing_data_proportion.csv")
data = data[(data$BBCH_50plus_proportion==1),]
str(data)

#create a column called budbreak_## to plot budbreak against time -- BE SURE YOUR COLUMN SUMMATIONS ARE CORRECT.
data=data%>%
  mutate(budbreak_51 = rowSums(data[,13:21]))

data=data%>%
  mutate(budbreak_52 = rowSums(data[,14:21]))

data = data %>%
  mutate(budbreak_53 = rowSums(data[,15:21]))

data=data%>%
  mutate(budbreak_54 = rowSums(data[,16:21]))
str(data)

#can look at budbreak (stage53+) at each date.
#Nov 16th, 2019 collection
ggplot(dplyr::filter(data, chill_portions_since_oct1=="27.8"), aes(x=hours_at_22.2C_field_adjusted, y=budbreak_53))+
  geom_point(aes(color=Genotype, shape=Genotype), size=4) +
  geom_line(aes(color=Genotype))+
  ggtitle("Nov 16th Collection") +
  scale_color_manual(values=c("lightpink", "lightpink", "lightpink", "steelblue", "steelblue", "steelblue", "steelblue", "steelblue"),limits = c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name = "Bloom Group")+
  scale_shape_manual(values=(c(19,19,19,15,15,15,15,15)), limits=c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name="Bloom Group")+
  xlab("Hours in heat (22.2C)") +
  ylab("Budbreak % (Stage 53+)") +
  ylim(0,1)+
  theme_minimal() +
  theme(plot.title=element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

#Dec 12th, 2019 collection
ggplot(dplyr::filter(data, chill_portions_since_oct1=="43.5"), aes(x=hours_at_22.2C_field_adjusted, y=budbreak_53))+
  geom_point(aes(color=Genotype, shape=Genotype), size=4) +
  geom_line(aes(color=Genotype))+
  ggtitle("Dec 12th Collection") +
  scale_color_manual(values=c("lightpink", "lightpink", "lightpink", "steelblue", "steelblue", "steelblue", "steelblue", "steelblue"),limits = c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name = "Bloom Group")+
  scale_shape_manual(values=(c(19,19,19,15,15,15,15,15)), limits=c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name="Bloom Group")+
  xlab("Hours in heat (22.2C)") +
  ylab("Budbreak % (Stage 53+)") +
  ylim(0,1)+
  theme_minimal() +
  theme(plot.title=element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

#Jan 10th, 2019 collection
ggplot(dplyr::filter(data, chill_portions_since_oct1=="60.6"), aes(x=hours_at_22.2C_field_adjusted, y=budbreak_53))+
  geom_point(aes(color=Genotype, shape=Genotype), size=4) +
  geom_line(aes(color=Genotype))+
  ggtitle("Jan 10th Collection") +
  scale_color_manual(values=c("lightpink", "lightpink", "lightpink", "steelblue", "steelblue", "steelblue", "steelblue", "steelblue"),limits = c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name = "Bloom Group")+
  scale_shape_manual(values=(c(19,19,19,15,15,15,15,15)), limits=c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name="Bloom Group")+
  xlab("Hours in heat (22.2C)") +
  ylab("Budbreak % (Stage 53+)") +
  ylim(0,1)+
  theme_minimal() +
  theme(plot.title=element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))
  

#Jan 28th, 2019 collection
ggplot(dplyr::filter(data, chill_portions_since_oct1=="64.8"), aes(x=hours_at_22.2C_field_adjusted, y=budbreak_53))+
  geom_point(aes(color=Genotype, shape=Genotype), size=4) +
  geom_line(aes(color=Genotype))+
  ggtitle("Jan 28th Collection") +
  scale_color_manual(values=c("lightpink", "lightpink", "lightpink", "steelblue", "steelblue", "steelblue", "steelblue", "steelblue"),limits = c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name = "Bloom Group")+
  scale_shape_manual(values=(c(19,19,19,15,15,15,15,15)), limits=c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name="Bloom Group")+
  xlab("Hours in heat (22.2C)") +
  ylab("Budbreak % (Stage 53+)") +
  ylim(0,1)+
  theme_minimal() +
  theme(plot.title=element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))


#Feb 6th, 2019 collection
ggplot(dplyr::filter(data, chill_portions_since_oct1=="69.3"), aes(x=hours_at_22.2C_field_adjusted, y=budbreak_53))+
  geom_point(aes(color=Genotype, shape=Genotype), size=4) +
  geom_line(aes(color=Genotype))+
  ggtitle("Feb 6th Collection") +
  scale_color_manual(values=c("lightpink", "lightpink", "lightpink", "steelblue", "steelblue", "steelblue", "steelblue", "steelblue"),limits = c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name = "Bloom Group")+
  scale_shape_manual(values=(c(19,19,19,15,15,15,15,15)), limits=c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name="Bloom Group")+
  xlab("Hours in heat (22.2C)") +
  ylab("Budbreak % (Stage 53+)") +
  ylim(0,1)+
  theme_minimal() +
  theme(plot.title=element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))


#Feb 17th, 2019 collection
ggplot(dplyr::filter(data, chill_portions_since_oct1=="69.9"), aes(x=hours_at_22.2C_field_adjusted, y=budbreak_53))+
  geom_point(aes(color=Genotype, shape=Genotype), size=4) +
  geom_line(aes(color=Genotype))+
  ggtitle("Feb 17th Collection") +
  scale_color_manual(values=c("lightpink", "lightpink", "lightpink", "steelblue", "steelblue", "steelblue", "steelblue", "steelblue"),limits = c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name = "Bloom Group")+
  scale_shape_manual(values=(c(19,19,19,15,15,15,15,15)), limits=c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name="Bloom Group")+
  xlab("Hours in heat (22.2C)") +
  ylab("Budbreak % (Stage 53+)") +
  ylim(0,1)+
  theme_minimal() +
  theme(plot.title=element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

#Feb 25th, 2019 collection
ggplot(dplyr::filter(data, chill_portions_since_oct1=="76.4"), aes(x=hours_at_22.2C_field_adjusted, y=budbreak_53))+
  geom_point(aes(color=Genotype, shape=Genotype), size=4) +
  geom_line(aes(color=Genotype))+
  ggtitle("Feb 25th Collection") +
  scale_color_manual(values=c("lightpink", "lightpink", "lightpink", "steelblue", "steelblue", "steelblue", "steelblue", "steelblue"),limits = c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name = "Bloom Group")+
  scale_shape_manual(values=(c(19,19,19,15,15,15,15,15)), limits=c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name="Bloom Group")+
  xlab("Hours in heat (22.2C)") +
  ylab("Budbreak % (Stage 53+)") +
  ylim(0,1)+
  theme_minimal() +
  theme(plot.title=element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

#Mar 1st, 2019 collection - these were kept in a cold room after collection and given 'artificial chilling'
ggplot(dplyr::filter(data, chill_portions_since_oct1=="91.5"), aes(x=hours_at_22.2C_field_adjusted, y=budbreak_53))+
  geom_point(aes(color=Genotype, shape=Genotype), size=4) +
  geom_line(aes(color=Genotype))+
  ggtitle("Mar 1st Collection, Artificially Chilled") +
  scale_color_manual(values=c("lightpink", "lightpink", "lightpink", "steelblue", "steelblue", "steelblue", "steelblue", "steelblue"),limits = c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name = "Bloom Group")+
  scale_shape_manual(values=(c(19,19,19,15,15,15,15,15)), limits=c('27-02-08'='27-02-08','27-03-08'='27-03-08','27-04-12'='27-04-12','27-02-65'='27-02-65','27-03-25'='27-03-25','27-03-28'='27-03-28','27-03-27'='27-03-27','27-03-46'='27-03-46'), name="Bloom Group")+
  xlab("Hours in heat (22.2C)") +
  ylab("Budbreak % (Stage 53+)") +
  ylim(0,1)+
  theme_minimal() +
  theme(plot.title=element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_text(face="bold", size=15), plot.title=element_text(size=18, face = "bold", hjust=0.5), legend.text=element_text(size=12), legend.title=element_text(size=15))

#################### you can do this for the other budbreak stages too :)


