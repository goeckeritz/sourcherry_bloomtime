################################
#                              #
# Learning how to use ChillR   #
#                              #
################################

setwd("/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/cherry_stuff_Charity/2019_winter_experiments/chillR/")

#First, download the temperature data and get it into the format described in 'chillR documentation.pdf'

#Install and load the package and its dependencies

#install.packages('chillR') 
#Everything properly downloaded - huzzah!


#Read in the data. I downloaded it from Clarksville's enviroweather. Temp is in C. Some editing was done in excel to get the date in that format, and the hours were changes in a text editer using regular expressions. 

data <- read.csv("oct1_2017_to_feb26_2022.csv", header=T)
str(data)
#Need to change the JDay to Julian Day of that year. It probably shouldn't be a factor either but we'll find out soon enough. 
#This package I found on Stack Overflow; let's see if I can get it to work and replace my current JDay column
#install.packages('lubridate')

require(lubridate)
require(chillR)
require(tidyverse)

data2 <- data %>%
  mutate(JDay = yday(JDay)) # you MUST CHANGE YOUR DATE TO YEAR-MON-DAY before this works, and have it open in excel, because it auto changes back the date format. God fuck excel. 
#Did it work?
head(data2)
str(data2)

###########################checking data format with sample data in chillR package################
#Yes! Now I need to make sure that each of the 4 variables 
#(columns) are the right kind of class.  
?chilling
?stack_hourly_temps
#I need to see an example data set that stack_hourly_temps produces
weather<-fix_weather(KA_weather[which(KA_weather$Year>2004),])
hourtemps<-stack_hourly_temps(weather, latitude=50.4)
str(hourtemps)
#Ah, so this is a list with two elements: hourtemps, and QC. Looks like Year is an integer, JDay is a num, Hour is num but I may have it in the wrong format initially (have 01 but should start with 0), and Temp is numeric. So all numeric except Year. 
data_example <- hourtemps$hourtemps
str(data_example)
range(data_example$Hour)
#indeed, the first hour is actually hour '0' and the last is '23'.

#simple to change!
############################################### 

data3 <- data2 %>%
  mutate(Hour = Hour - 1)
str(data3)
range(data3$Hour)
range(data3$Temp) #idk why R is angry about this, chillR calculates chill just fine. 
head(data3)
tail(data3)
#We all set, boiiiiiiii

calculations <- chilling(THourly = data3, Start_JDay = 1, End_JDay = 365)
#It worked! Woo! Super easy and super useful!!!! Utah model calculations are odd, but there's no limit on the negations for high temp so that's likely why. 
#What if I wanted to know how many chill portions I had between Oct 1st and Oct 13th, which is when my first collection date was?
oct1_to_oct13 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 287) # 6.6 chill portions
#chill since my second collection, all years?
oct1_to_nov16 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 320) #27.8 chill portions
#chill since my most recent (third) collection, all years?
oct1_to_dec12 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 347) # 43.5 chill portions
oct1_to_dec17 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 352)
oct1_to_dec31 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 365) # 53.5 chill portions
oct1_to_jan10 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 10) # 60.6 chill portions
oct1_to_jan28 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 28) # 64.8 chill portions
oct1_to_feb6 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 37) #69.3 chill portions
oct1_to_feb17 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 48) #69.9 chill portions
oct1_to_feb25 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 56) # 76.4 chill portions
oct1_to_mar1  <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 61) #78.9 chill portions 
oct1_to_mar17  <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 77) #artifically induced chill from mar1 - mar17 # 91.5 chill portions 
oct1_to_oct29 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 302) #8.9 chill portions for 2021
oct1_to_nov14 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 319) #20.7 cp for 2021
oct1_to_nov24 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 329) #28 cp for 2021
oct1_to_dec4 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 339)
oct1_to_dec11 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 346)
oct1_to_dec29 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 364) #46.3 cp for 2021
oct1_to_jan5 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 5) #
oct1_to_jan13 <- chilling(THourly = data3, Start_JDay = 274, End_JDay = 13) #50 cp for 2021-22
oct1_to_jan27 <- chilling(THourly = data3, Start_JDay = 274, End_JDay=27) #54 cp for 2021-2022 
oct1_to_feb10 <- chilling(THourly = data3, Start_JDay = 274, End_JDay=41) #60.8 cp for 2021-2022
oct1_to_feb26 <- chilling(THourly = data3, Start_JDay = 274, End_JDay=57)  #67.7 cp for 2021-2022
jan1_to_feb22 <- chilling(THourly = data3, Start_JDay = 1, End_JDay=53)


#Negligible chill portions since Dec12. There's likely NO need to collect chill branches this Friday!




#######################################################################
#these calculations for chill portions are for the ROS collection dates. 
#All temp data is from natural field conditions. Note that a gap in hourly data near 2/14/2020
#led me to substitute those values by the next nearest weather station. 

setwd("C://Users/Goeckeritz/Desktop/cherry_stuff_Charity/2019_winter_experiments/chillR/")

require(lubridate)
require(chillR)
require(tidyverse)

data5 <- read.csv("oct1_2017_to_april1_2020.csv", header=T)
str(data5)

data5 = data5 %>%
  mutate(JDay = yday(JDay)) %>%
  mutate(Hour = Hour - 1)

str(data5) #all set. 

oct1_to_oct13 <- chilling(THourly = data5, Start_JDay = 274, End_JDay = 287) # 6.6 chill portions
oct1_to_nov16 <- chilling(THourly = data5, Start_JDay = 274, End_JDay = 320) # 26.6 chill portions
oct1_to_dec20 <- chilling(THourly = data5, Start_JDay = 274, End_JDay = 354) # 43.7 chill portions
oct1_to_feb22 <- chilling(THourly = data5, Start_JDay = 274, End_JDay = 53) # 74.8 chill portions
oct1_to_mar23 <- chilling(THourly = data5, Start_JDay = 274, End_JDay = 83) #  chill portions #JDay 83 cuz 2020 was a leap year. 



