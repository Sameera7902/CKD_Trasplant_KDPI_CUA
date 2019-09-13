# HE_Preprocessing waitlist data and dialysis data
# process the data ready for analysis - Markov model Scenario 02
# August 2019



# Load libraries
library(foreign)
library(dplyr)
library(readr)
library(lubridate)
library(tidyverse)
library(data.table)
library(forcats)
library(survival)
library(eha)
library(SurvRegCensCov) # Intepret Weibull parameters
library(survminer)
library(ggplot2)

#### Import the dialysis data base
setwd("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA")

load("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA/Dialysis_processedData.RData")
waitlist_dialysis = read.csv("Waitlist_dialysis.csv")
transpData_dialysis = read.csv("TranspIN_DialysisDara.csv")



### Parse as dates
waitlist_dialysis$waitdate = dmy(waitlist_dialysis$waitdate)

### Select the important columns
waitlist_dialysis = waitlist_dialysis %>%
  select(1, 3,4,6)


# Joining data sets
waitlist_dialysis = data.table(waitlist_dialysis)
transpData_dialysis = data.table(transpData_dialysis)
dialysis = data.table(dialysis) 

dial_trans = merge(dialysis, transpData_dialysis, by = "id", all.x = TRUE ) #Combine dialysis and transplant data sets 



dial_trans = dial_trans %>%
  select(1:10)


# excluding transplanted patients
dial_NOtran = dial_trans %>%
  filter(is.na(graftno))

dial_NOtran = dial_NOtran %>%
  select(1:9)


# Joining the waitlist data
dialy_merged = merge(dial_NOtran, waitlist_dialysis, by = "id", allow.cartesian = TRUE)

# No of waitlist in the "waitlist_trans+dilaysis" data set
waitlist_dialysis%>%
  filter(waitseq == "1") %>%
  summarise(n())

# No of waitlist in the "waitlsted_no transplant_dialysis" data set
dialy_merged%>%
  filter(waitseq == "1") %>%
  summarise(n())

dialy_merged_1 = dialy_merged %>%
  filter(waitseq == "1")

# For PSA
dialy_merged_1 %>%
  filter(waitseq == "1" & ageatRRTstart >39 & ageatRRTstart < 50) %>% # KDPI < 25
  summarise(prop = (mean(death == "Death")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


dialy_merged_1 %>%
  filter(waitseq == "1" & ageatRRTstart > 49 & ageatRRTstart < 60) %>% # KDPI < 25
  summarise(prop = (mean(death == "Death")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)

dialy_merged_1 %>%
  filter(waitseq == "1" & ageatRRTstart > 59 & ageatRRTstart < 70) %>% # KDPI < 25
  summarise(prop = (mean(death == "Death")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)




# Kaplan-Meier plot
plot_M = plot(survival::survfit(Surv(duration, death == "Death") ~ 1, data =dialy_merged_1 ), 
                      main = "Survival while waitlisted for transplant",
                      xlab = "Years", 
                      ylab = "Overall survival probability",
                      )


wbmod_waitlist = survreg(Surv(duration, death == "Death") ~ ageatRRTstart, data = dialy_merged_1, dist='weibull') # Gender was not significant in the model.  Thus removed
summary(wbmod_waitlist)
ConvertWeibull(wbmod_waitlist,conf.level = 0.95)







