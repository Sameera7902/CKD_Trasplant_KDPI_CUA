# Probability calculation
# process the data ready for analysis - Markov model Scenario 01 (Dialysis data set)
# July 2019

# Load libraries
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)
library(survival)
library(SurvRegCensCov) # Intepret Weibull parameters
library(survminer)
library(ggplot2)
library(eha)

#### Import the dialysis data base
setwd("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA")
load("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA/Dialysis_Transplant_merged_ReadyforAnalysis.RData")

# Kaplan-Meier plot
plot_dialysis_mortality = plot(survival::survfit(Surv(duration, death == "Death") ~ 1, data =dialysis_transplant_merged), 
                      main = "Dialysis - Mortality",
                      xlab = "Years", 
                      ylab = "Overall survival probability",
                     )


### Weibull Regresson
dialysis_transplant_merged = mutate(dialysis_transplant_merged, duration = duration + 0.001) # The survreg function in R does not allow time = 0


wbmod_dialysis_M = survreg(Surv(duration, death == "Death") ~ ageatRRTstart, data = dialysis_transplant_merged, dist='weibull') 
summary(wbmod_dialysis_M )
ConvertWeibull(wbmod_dialysis_M ,conf.level = 0.95)

#$$ Testing Weibull vs Parametric model
wbmod_dialysis_M_1 = phreg(Surv(duration, death == "Death") ~ ageatRRTstart, data = dialysis_transplant_merged, dist='weibull')
dialysis_M_Cox = coxreg(Surv(duration, death == "Death") ~ ageatRRTstart, data = dialysis_transplant_merged)
check.dist(wbmod_dialysis_M_1 , dialysis_M_Cox)

AIC(wbmod_dialysis_M) # Smaller the value, better the model
AIC(dialysis_M_Cox)

BIC(wbmod_dialysis_M)
BIC(dialysis_M_Cox)


### Exponential Regresson
exmod_dialysis_M = survreg(Surv(duration, death == "Death") ~ ageatRRTstart, data = dialysis_transplant_merged, dist='exponential') 
summary(exmod_dialysis_M ) #? Lambda value



