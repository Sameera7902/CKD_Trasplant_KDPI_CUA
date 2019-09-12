# HE_Preprocessing
# process the data ready for analysis - Markov model Scenario 01 (Dialysis data set)
# July 2019


# Load libraries
#install.packages("foreign")
library(foreign)
library(dplyr)
library(readr)
#install.packages("lubridate")
library(lubridate)
#install.packages("tidyr")
library(tidyr)

#### Import the dialysis data base
setwd("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA")
load("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA/Dialysis_Transplant_merged.RData")

#Remove transplant patient data
dialysis_transplant_merged = dialysis_transplant_merged %>%
  filter(is.na(ageattransplant))

dialysis_transplant_merged = dialysis_transplant_merged %>%
  select(1:9)


# Create death column (1 = Death; 0 = Alive)
dialysis_transplant_merged = dialysis_transplant_merged %>%
  mutate(death = ifelse(ageatdeath.x >0, "Death", "Alive"))

dialysis_transplant_merged$death = ifelse(is.na(dialysis_transplant_merged$death), "Alive", dialysis_transplant_merged$death)

dialysis_transplant_merged %>%
  summarise(prop = mean(death == "Death"))

dialysis_transplant_merged$death = as.factor(dialysis_transplant_merged$death)


# Varibale convertion
dialysis_transplant_merged$deathdate.x = as.character(dialysis_transplant_merged$deathdate.x)

# Parse as dates
dialysis_transplant_merged$deathdate.x = parse_date(dialysis_transplant_merged$deathdate.x, format = "%d %b %Y")


# Create duration of dialysis column
dialysis_transplant_merged = dialysis_transplant_merged %>%
  mutate(duration =  deathdate.x - treatmentdate.x)

dialysis_transplant_merged$duration = ifelse(is.na(dialysis_transplant_merged$duration), (as.Date("2017-12-31")- dialysis_transplant_merged$treatmentdate.x), dialysis_transplant_merged$duration)

dialysis_transplant_merged = dialysis_transplant_merged %>%
  mutate(duration = duration/365)

# save the data
save(dialysis_transplant_merged, file='Dialysis_Transplant_merged_ReadyforAnalysis.RData')


