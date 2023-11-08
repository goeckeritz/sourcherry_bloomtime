#########Plotting the sectioning staging for 2018-2019, and 2019-2020

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/Pop4_bloomtime_PAPER/stages/")

library(ggplot2)
library(tidyverse)
library(dplyr)
library(tibble)

#the dataaaaa
#real simple -- i wanna make a plot with both years on it, as well as separate years. 
#do year as shape, and different shades of blue for earlies and different shades of pink for lates. 
#update 11/1/2023 - Rebecca recommended doing a Mann-Whitney test fo our discrete staging data. Data can't be
#normal since it's discrete. We do have an ordinal measurement scale too. https://www.statisticssolutions.com/free-resources/directory-of-statistical-analyses/mann-whitney-u-test/ 

data <- read.table("staging_cherry_data.txt", header=TRUE, stringsAsFactors = FALSE)
str(data)

data <- data %>%
  select(year, julian_day, Julian_days_since_Jun_1st, bloom_group, individual, stage)

str(data)

data$year = as.factor(data$year)
data$bloom_group = as.factor(data$bloom_group)
data$individual = as.factor(data$individual)
data$stage = ordered(as.factor(data$stage))
levels(data$stage)
is.ordered(data$stage)
data$stage_numeric = as.numeric(as.character(data$stage))
data$tree_id = as.factor(paste(data$bloom_group, data$individual, sep="_"))

str(data)
levels(data$tree_id)

#both years
ggplot(data, aes(Julian_days_since_Jun_1st, stage, color=tree_id, shape=year)) +
  geom_point(position=position_jitter(w=1, h=0.5),size=3)+
  ggtitle("Development Based on Histological Sections") +
  xlab("Days since June 1st") +
  ylab("Developmental Stage") + 
  scale_color_manual(values=c("steelblue", "steelblue2", "steelblue1", "skyblue", "lightblue2", "pink2", "violetred2", "violetred3", "violetred4")) +
  scale_shape_manual(values=c(1,16)) +
  theme_minimal()

ggplot(data=data, aes(Julian_days_since_Jun_1st, stage, color=tree_id)+
               geom_point(position=position_jitter(w=1, h=0.5),size=3))+
               ggtitle("Development Based on Histological Sections") +
               xlab("Days since June 1st") +
               ylab("Developmental Stage") + 
               scale_color_manual(values=c("steelblue", "steelblue2", "steelblue1", "skyblue", "lightblue2", "pink2", "violetred2", "violetred3", "violetred4")) +
               scale_shape_manual(values=c(1,16)) +
               theme_minimal()

data_2018 = data %>%
  dplyr::filter(data$year!="2019-20")

data_2019 = data %>%
  dplyr::filter(data$year!="2018-19")

ggplot(data_2018, aes(Julian_days_since_Jun_1st, stage, color=tree_id, shape=tree_id))+
  geom_point(position=position_jitter(w=2, h=0.2),size=4.5, stroke=1)+
  ggtitle("Development Based on Histological Sections (2018-19)") +
  xlab("Days since June 1st, 2018") +
  ylab("Developmental Stage") + 
  scale_color_manual(values=c("steelblue", "steelblue2", "steelblue1", "pink2", "violetred2", "violetred3"), labels=c("early_1" = "Early 1", "early_2"="Early 2", "early_3"="Early 3", "late_1" = "Late 1", "late_2" = "Late 2", "late_3" = "Late 3"), name = "Bloom Group") +
  scale_shape_manual(values=c(0, 0, 0, 1, 1, 1), labels=c("early_1" = "Early 1", "early_2"="Early 2", "early_3"="Early 3", "late_1" = "Late 1", "late_2" = "Late 2", "late_3" = "Late 3"), name = "Bloom Group") +
  scale_x_continuous(limits=c(0, 350), breaks=c(0,30,61,92,122,153,183,214,245,273,304,334),labels=c("Jun", "Jul", "Aug", "Sep", "Oct","Nov","Dec", "Jan", "Feb", "Mar", "Apr","May")) +
  scale_y_discrete(limits=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"), breaks=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"), labels=c("","","2","","4","","6","","8","","10","","12","","14","","16","","18","","20",""))+
  theme_minimal()+
  theme(axis.title.x=element_text(face="bold", size=18)) + theme(axis.title.y=element_text(face="bold", size=18)) + theme(legend.title = element_text(face="bold", size=18)) + theme(legend.text = element_text(size=15, face="bold")) + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 15, face="bold"), axis.text.y=element_text(size=15, face="bold"), title=element_text(size=15, face="bold"))


ggplot(data_2019, aes(Julian_days_since_Jun_1st, stage, color=tree_id, shape=tree_id))+
  geom_point(position=position_jitter(w=2, h=0.2),size=4.5, stroke=1)+
  ggtitle("Development Based on Histological Sections (2019-20)") +
  xlab("Days since June 1st, 2019") +
  ylab("Developmental Stage") + 
  scale_color_manual(values=c("steelblue", "steelblue2", "steelblue1", "skyblue", "lightblue2", "pink2", "violetred2", "violetred3", "violetred4"), labels=c("early_1" = "Early 1", "early_2"="Early 2", "early_3"="Early 3", "early_4"= "Early 4", "early_5" = "Early 5", "late_1" = "Late 1", "late_2" = "Late 2", "late_3" = "Late 3", "late_4" = "Late 4"), name = "Bloom Group") +
  scale_shape_manual(values=c(0, 0, 0, 0, 0, 1, 1, 1, 1), labels=c("early_1" = "Early 1", "early_2"="Early 2", "early_3"="Early 3", "early_4"= "Early 4", "early_5" = "Early 5", "late_1" = "Late 1", "late_2" = "Late 2", "late_3" = "Late 3", "late_4" = "Late 4"), name = "Bloom Group") +
  scale_x_continuous(limits=c(0, 350), breaks=c(0,30,61,92,122,153,183,214,245,273,304,334),labels=c("Jun", "Jul", "Aug", "Sep", "Oct","Nov","Dec", "Jan", "Feb", "Mar", "Apr","May")) +
  scale_y_discrete(limits=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"), breaks=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"), labels=c("","","2","","4","","6","","8","","10","","12","","14","","16","","18","","20",""))+
  theme_minimal()+
  theme(axis.title.x=element_text(face="bold", size=18)) + theme(axis.title.y=element_text(face="bold", size=18)) + theme(legend.title = element_text(face="bold", size=18)) + theme(legend.text = element_text(size=15, face="bold")) + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 15, face="bold"), axis.text.y=element_text(size=15, face="bold"), title=element_text(size=15, face="bold"))


##Mann-Whitney tests by date...

#2018-2019:
#Goal: at pre-dormancy dates, do early stages and late stages come from the same distribution? I.e.,
#are they statistically different from one another? Not sure I can do Nesting here. Found this: https://rdrr.io/cran/nestedRanksTest/man/nestedRanksTest.html
#But I might blow the code up if I use an ordinal dependent. The description seems to fit the data.

#we'll test dates 220, 263, and 303.(August 8, Sep 20, Oct 30) --- Excluding Jun14 cuz... duh, they're all at the same stage. 

Aug8_2018 = data_2018 %>%
  dplyr::filter(data_2018$julian_day==220)

Sept20_2018 = data_2018 %>%
  dplyr::filter(data_2018$julian_day==263)

Oct30_2018 = data_2018 %>%
  dplyr::filter(data_2018$julian_day==303)

Nov28_2018 = data_2018 %>%
  dplyr::filter(data_2018$julian_day==332)

###

Jul18_2019 = data_2019 %>%
  dplyr::filter(data_2019$julian_day==199)

Aug8_2019 = data_2019 %>%
  dplyr::filter(data_2019$julian_day==220)

Aug22_2019 = data_2019 %>%
  dplyr::filter(data_2019$julian_day==234)

Sep7_2019 = data_2019 %>%
  dplyr::filter(data_2019$julian_day==250)

Nov16_2019 = data_2019 %>%
  dplyr::filter(data_2019$julian_day==320)

#After some discussion with Andrea, we decided -- although it's not optimal -- that we'd treat each flower bud as an independent observation for 
#each bloom group, because there isn't an easy way with the Mann-Whitney test to control for the random effect of each tree. 
#The only other option I could come up with was test the medians of each tree for each date... but that would completely destroy my
#statistical power, to the point that it's likely not worth to do. 

#so, we'll do it as described and consider a bonferroni correction for multiple tests. Let's start with 2018-19
#We've 4 dates to test: 
#Aug8, Sep20, and Oct30

#Then we'll do 2019-2020. We have X dates to test:
#Jul18, Aug8, Aug22, Sep9, Nov13

#Question... do I have to coerce my ordinal data into an integer? Is that... okay??
#"While both are nonparametric and involve summation of ranks, the Wilcoxon signed-rank test requires that the data is paired while the Wilcoxon rank-sum test is used for unpaired data."
#https://en.wikipedia.org/wiki/Talk:Wilcoxon_signed-rank_test#Ordinal_data

#this example seems pretty similar to mine, and these guys just willy nilly start treating the data as numeric
#with no EFFING EXPLANATION: https://rcompanion.org/handbook/F_04.html
#Lord of all statistics, WHY???
#In general I can say that there is some numerical value to the stagings, sure -- 1 < 2 < 3.... etc.
#But I'm definitely not confident about saying the difference between 1 and 2 is the same as the difference between 2 and 3. 
#I'll treat it as numeric and consult the co-authors. 

install.packages("coin")

str(Aug8_2018)

##For 2018-2019 pre-dormancy season: 

wilcox.test(stage_numeric ~ bloom_group, data=Aug8_2018) ## uncorrected p-value: 0.02747 ||| Bonferroni corrected p-value: 0.247
wilcox.test(stage_numeric ~ bloom_group, data=Sept20_2018) ## uncorrected p-value: 0.007526 ||| Bonferroni corrected p-value: 0.0677
wilcox.test(stage_numeric ~ bloom_group, data=Oct30_2018) ## uncorrected p-value: 0.008678 ||| Bonferroni corrected p-value: 0.0781
wilcox.test(stage_numeric ~ bloom_group, data=Nov28_2018) ## uncorrected p-value: 0.03015 ||| Bonferroni corrected p-value: 0.271

##For 2019-2020 pre-dormancy season:

wilcox.test(stage_numeric ~ bloom_group, data=Jul18_2019) ## uncorrected p-value: 0.09426 ||| Bonferroni corrected p-value: 0.84
wilcox.test(stage_numeric ~ bloom_group, data=Aug8_2019) ## uncorrected p-value: 0.005488 ||| Bonferroni corrected p-value: 0.0494
wilcox.test(stage_numeric ~ bloom_group, data=Aug22_2019) ## uncorrected p-value: 0.00001996 ||| Bonferroni corrected p-value: 0.00018
wilcox.test(stage_numeric ~ bloom_group, data=Sep7_2019) ## uncorrected p-value: 0.00007961 ||| Bonferroni corrected p-value: 0.000716
wilcox.test(stage_numeric ~ bloom_group, data=Nov16_2019) ## uncorrected p-value: 0.1138 ||| Bonferroni corrected p-value: 1


#There are 9 total tests.. so a classical Bonferroni would be to simply multiply all p-values by 9. 
#unsure if this will make it into the manuscript yet. 














             




