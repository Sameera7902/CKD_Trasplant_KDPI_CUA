# Rcode_HE_Sc1_transpalnt
# Calculate the probabilities for the Markov model - Probabilistic Sensitivity Analysis
# July 2019


# get the data
# Office Comp
setwd("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA")
load("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA/HE_Sc1_transplant.RData")
load("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA/Dialysis_Transplant_merged_ReadyforAnalysis.RData")



# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(survival)
library(eha)
library(SurvRegCensCov) # Intepret Weibull parameters
library(survminer)
library(ggplot2)
#install.packages("gmodels")
library(gmodels)

### Pre-processing
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, surv_years_1 = surv_years + 0.001) # The survreg function in R does not allow time = 0
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, graft_years_1 = graft_years + 0.001) # The survreg function in R does not allow time = 0


# Distribution of "age at transplant"
summarise(transplant_deceased_cl_selected$ageattransplant)




# Creating separate data sets
kdpi_25 = filter(transplant_deceased_cl_selected, kdpi < 25 & graftno == "1")
kdpi_50 = filter(transplant_deceased_cl_selected, kdpi >= 25 & kdpi < 50 & graftno == "1")
kdpi_75 = filter(transplant_deceased_cl_selected, kdpi  >= 50 & kdpi < 75 & graftno == "1")
kdpi_100 = filter(transplant_deceased_cl_selected, kdpi >= 75 & graftno == "1")



# @@@@@@@@@ Mortality following the 1st transplant - 40 year analysys ( Last 5 yr data) @@@@@@@@@

# KDPI <25
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 39 & ageattransplant < 50 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 25-49
transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 39 & ageattransplant < 50 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)

# KDPI 50-74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 39 & ageattransplant < 50 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)

# KDPI > 74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 39 & ageattransplant < 50 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Mortality following the 1st transplant - 40 year analysys ( Total data) @@@@@@@@@

# KDPI <25
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 39 & ageattransplant < 50 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 25-49
transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 39 & ageattransplant < 50 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 50-74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 39 & ageattransplant < 50 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI > 74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 39 & ageattransplant < 50 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Graft falure following the 1st transplant - 40 year analysis (Last 5 yr data) @@@@@@@@@
# KDPI <25
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 39 & ageattransplant < 50 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 25-49
transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 39 & ageattransplant < 50 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)

# KDPI 50-74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 39 & ageattransplant < 50 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)

# KDPI > 74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 39 & ageattransplant < 50 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Graft falure following the 1st transplant - 40 year analysis (Total data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 39 & ageattransplant < 50) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 39 & ageattransplant < 50) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 39 & ageattransplant < 50 ) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 39 & ageattransplant < 50 ) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)




# @@@@@@@@@ Mortality following the re transplant - 40 years ( 5 year data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(graftno != "1" & ageattransplant > 39 & ageattransplant < 50 & transplantdate > "2012-01-01") %>% 
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# @@@@@@@@@ Mortality following the re transplant - 40 years ( All data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(graftno != "1" & ageattransplant > 39 & ageattransplant < 50 ) %>% 
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Mortality following dialysis - 40 year analysis @@@@@@@@@
dialysis_transplant_merged %>%
  filter(ageatRRTstart > 39 & ageatRRTstart < 50) %>% 
  summarise(prop = (mean(death == "Death")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)






##################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$









# @@@@@@@@@ Mortality following the 1st transplant - 50 year analysys ( Last 5 yr data) @@@@@@@@@

# KDPI <25
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 49 & ageattransplant < 60 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 25-49
transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 49 & ageattransplant < 60 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 50-74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 49 & ageattransplant < 60 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI > 74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 49 & ageattransplant < 60 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# @@@@@@@@@ Mortality following the 1st transplant - 50 year analysys ( Total data) @@@@@@@@@

# KDPI <25
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 49 & ageattransplant < 60 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 25-49
transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 49 & ageattransplant < 60 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 50-74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 49 & ageattransplant < 60 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI > 74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 49 & ageattransplant < 60 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Graft falure following the 1st transplant - 50 year analysis (Last 5 yr data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 49 & ageattransplant < 60 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 49 & ageattransplant < 60 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 49 & ageattransplant < 60 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 49 & ageattransplant < 60 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Graft falure following the 1st transplant - 50 year analysis (Total data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 49 & ageattransplant < 60) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 49 & ageattransplant < 60 ) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 49 & ageattransplant < 60 ) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 49 & ageattransplant < 60 ) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)




# @@@@@@@@@ Mortality following the re transplant - 50 years ( 5 year data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(graftno != "1" & ageattransplant > 49 & ageattransplant < 60 & transplantdate > "2012-01-01") %>% 
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# @@@@@@@@@ Mortality following the re transplant - 50 years ( All data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(graftno != "1" & ageattransplant > 49 & ageattransplant < 60 ) %>% 
  summarise(prop = (mean(lastknownstatus == "Dead")),
           Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Mortality following dialysis - 50 year analysis @@@@@@@@@
dialysis_transplant_merged %>%
  filter(ageatRRTstart > 49 & ageatRRTstart < 60) %>% 
  summarise(prop = (mean(death == "Death")),
           Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)








# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$








# @@@@@@@@@ Mortality following the 1st transplant - 60 year analysys ( Last 5 yr data) @@@@@@@@@

# KDPI <25
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 59 & ageattransplant < 71 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# KDPI 25-49
transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 59 & ageattransplant < 71 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 50-74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 59 & ageattransplant < 71 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI > 74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & aageattransplant > 59 & ageattransplant < 71 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Mortality following the 1st transplant - 60 year analysys ( Total data) @@@@@@@@@

# KDPI <25
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 59 & ageattransplant < 71) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# KDPI 25-49
transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 59 & ageattransplant < 71 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)

# KDPI 50-74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 59 & ageattransplant < 71 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)

# KDPI > 74
transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 59 & ageattransplant < 71 ) %>% # KDPI < 25
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Graft falure following the 1st transplant - 60 year analysis (Last 5 yr data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 59 & ageattransplant < 71 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 59 & ageattransplant < 71 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 59 & ageattransplant < 71 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 59 & ageattransplant < 71 & transplantdate > "2012-01-01") %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
             Total = n(),
             z = qnorm(0.975),
             rate = (-log(1-prop))/5,
             prob = (1- exp(-rate)),
             se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Graft falure following the 1st transplant - 60 year analysis (Total data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & ageattransplant > 59 & ageattransplant < 71) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



transplant_deceased_cl_selected %>%
  filter(kdpi >= 25 & kdpi < 50 & ageattransplant > 59 & ageattransplant < 71 ) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 50 & kdpi < 75 & ageattransplant > 59 & ageattransplant < 71 ) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


transplant_deceased_cl_selected %>%
  filter(kdpi >= 75 & ageattransplant > 59  ) %>% # KDPI < 25
  summarise(prop = (mean(endtransplantcode == "Graft failure")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)




# @@@@@@@@@ Mortality following the re transplant - 60 years ( 5 year data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(graftno != "1" & ageattransplant > 59 & transplantdate > "2012-01-01") %>% 
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/5,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)


# @@@@@@@@@ Mortality following the re transplant - 60 years ( All data) @@@@@@@@@
transplant_deceased_cl_selected %>%
  filter(graftno != "1" & ageattransplant > 59 ) %>% 
  summarise(prop = (mean(lastknownstatus == "Dead")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# @@@@@@@@@ Mortality following dialysis - 60 year analysis @@@@@@@@@
dialysis_transplant_merged %>%
  filter(ageatRRTstart > 59 & ageatRRTstart < 71) %>% 
  summarise(prop = (mean(death == "Death")),
            Total = n(),
            z = qnorm(0.975),
            rate = (-log(1-prop))/10,
            prob = (1- exp(-rate)),
            se = sqrt(prob*(1-prob)/Total)) # standard error)



# KDPI vs Age at trasnsplant

plot_1 = ggplot(transplant_deceased_cl_selected, aes(x = (ageattransplant)))
plot_1 = plot_1 + geom_point(aes(y = kdpi, colour = "Kidney Transplants"))
plot_1 <- plot_1 + scale_colour_manual(values = c("orange"))
plot_1 <- plot_1 + labs(y = "KDPI value",
                        x = "Age at transplant",
                        colour = "",
                        title = "Age at transplant vs KDPI")
plot_1 <- plot_1 + theme(legend.position = c(0.125, 1))
plot_1 <- plot_1 + theme(plot.title = element_text(hjust = 0.5))
plot_1

plot_2 <- plot_1 + geom_hline(yintercept=25, color = "red", size=2)
plot_2

plot_3 <- plot_1 + geom_hline(yintercept=75, color = "red", size=2)
plot_3


            

