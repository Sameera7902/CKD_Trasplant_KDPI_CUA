# Rcode_HE_Sc1_transpalnt
# Calculate the probabilities for the Markov model using Weibull Regression
# July 2019


# get the data
# Office Comp
setwd("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA")
load("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/Scenario 1/R Work/HE_Sc1_transplant.RData")



# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(survival)
#install.packages("eha")
library(eha)
#install.packages("SurvRegCensCov")
library(SurvRegCensCov) # Intepret Weibull parameters
#install.packages("survminer")
library(survminer)
library(ggplot2)


#### @@@@@@@@@@@@ Data Preprocessing @@@@@@@@@@@@@@@@@

### Creating Survival years column
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, surv_years_1 = surv_years + 0.001) # The survreg function in R does not allow time = 0
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, graft_years_1 = graft_years + 0.001) # The survreg function in R does not allow time = 0


# Creating separate data sets based on KDPI levels
kdpi_25 = filter(transplant_deceased_cl_selected, kdpi < 25 & graftno == "1")
kdpi_50 = filter(transplant_deceased_cl_selected, kdpi >= 25 & kdpi < 50 & graftno == "1")
kdpi_75 = filter(transplant_deceased_cl_selected, kdpi  >= 50 & kdpi < 75 & graftno == "1")
kdpi_100 = filter(transplant_deceased_cl_selected, kdpi >= 75 & graftno == "1")


# Kaplan-Meier plot for different KDPI levels
km_kdpi_25_M = survfit(Surv(surv_years_1, lastknownstatus == "Dead") ~ 1, data = kdpi_25)
km_kdpi_50_M = survfit(Surv(surv_years_1, lastknownstatus == "Dead") ~ 1, data = kdpi_50)
km_kdpi_75_M = survfit(Surv(surv_years_1, lastknownstatus == "Dead") ~ 1, data = kdpi_75)
km_kdpi_100_M = survfit(Surv(surv_years_1, lastknownstatus == "Dead") ~ 1, data = kdpi_100)

plot_kdpi_25_M = plot(survival::survfit(Surv(surv_years_1, lastknownstatus == "Dead") ~ 1, data = kdpi_25), 
    main = "KDPI < 25",
    xlab = "Years", 
    ylab = "Overall survival probability",
    ylim=c(0.6,1))

plot_kdpi_50_M = plot(survival::survfit(Surv(surv_years_1, lastknownstatus == "Dead") ~ 1, data = kdpi_50), 
                    main = "KDPI 25-49",
                    xlab = "Years", 
                    ylab = "Overall survival probability",
                    ylim=c(0.6,1))

plot_kdpi_75_M = plot(survival::survfit(Surv(surv_years_1, lastknownstatus == "Dead") ~ 1, data = kdpi_75), 
                    xlab = "Years", 
                    ylab = "Overall survival probability",
                    ylim=c(0.6,1),
                    main = "KDPI 50-74")

plot_kdpi_100_M = plot(survival::survfit(Surv(surv_years_1, lastknownstatus == "Dead") ~ 1, data = kdpi_100), 
                    xlab = "Years", 
                    ylab = "Overall survival probability",
                    ylim=c(0.6,1),
                    main = "KDPI >74")



#### @@@@@@@@@@@@ Weibull Regresson - Transition probabilities @@@@@@@@@@@@@@@@@

#@@@@@@@ KDPI <25
wbmod_mortality_25 = survreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_25, dist='weibull') # Gender was not significant in the model.  Thus removed
summary(wbmod_mortality_25)
ConvertWeibull(wbmod_mortality_25,conf.level = 0.95)

#$$ Testing Weibull vs Parametric model
kdpi_25_WB = phreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_25, dist='weibull')
kdpi_25_Cox = coxreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_25)
check.dist(kdpi_25_WB,kdpi_25_Cox)

AIC(kdpi_25_Cox) # Smaller the value, better the model
AIC(wbmod_mortality_25)

BIC(kdpi_25_Cox)
BIC(wbmod_mortality_25)


#@@@@@@@ KDPI 25-49
wbmod_mortality_50 = survreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_50, dist='weibull') # Gender was not significant in the model.  Thus removed
summary(wbmod_mortality_50)
ConvertWeibull(wbmod_mortality_50,conf.level = 0.95)

#$$ Testing Weibull vs Parametric model
kdpi_50_WB = phreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_50, dist='weibull')
kdpi_50_Cox = coxreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_50)
check.dist(kdpi_50_WB,kdpi_50_Cox)

AIC(kdpi_50_Cox) # Smaller the value, better the model
AIC(wbmod_mortality_50)

BIC(kdpi_50_Cox)
BIC(wbmod_mortality_50)


#@@@@@@@ KDPI 50-74
wbmod_mortality_75 = survreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_75, dist='weibull') # Gender was not significant in the model.  Thus removed
summary(wbmod_mortality_75)
ConvertWeibull(wbmod_mortality_75,conf.level = 0.95)

#$$ Testing Weibull vs Parametric model
kdpi_75_WB = phreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_75, dist='weibull')
kdpi_75_Cox = coxreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_75)
check.dist(kdpi_75_WB,kdpi_75_Cox)

AIC(kdpi_75_Cox) # Smaller the value, better the model
AIC(wbmod_mortality_75)

BIC(kdpi_75_Cox)
BIC(wbmod_mortality_75)


#@@@@@@@ KDPI 75-100
wbmod_mortality_100 = survreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_100, dist='weibull') # Gender was not significant in the model.  Thus removed
summary(wbmod_mortality_100)
ConvertWeibull(wbmod_mortality_100,conf.level = 0.95)

#$$ Testing Weibull vs Parametric model
kdpi_100_WB = phreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_100, dist='weibull')
kdpi_100_Cox = coxreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = kdpi_100)
check.dist(kdpi_100_WB,kdpi_100_Cox)

AIC(kdpi_100_Cox) # Smaller the value, better the model
AIC(wbmod_mortality_100)

BIC(kdpi_100_Cox)
BIC(wbmod_mortality_100)


### Normal probabilities - not to be used in the mddel.  Just to check!

# Graft years = Less than 1 year
deaths_1 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat == "Less than 1 year" & lastknownstatus == "Dead")%>%
  summarise (n())
deaths_1

denom_1 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1")%>%
  summarise (n())
denom_1
  
p1 = deaths_1/denom_1
  
# Graft years = 1 - 2 years
deaths_2 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat == "1 - 2 years" & lastknownstatus == "Dead")%>%
  summarise (n())
deaths_2

denom_2 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat != "Less than 1 year")%>%
  summarise (n())
denom_2

p2 = deaths_2/denom_2

# Graft years = 2 - 3 years
deaths_3 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat == "2 - 3 years" & lastknownstatus == "Dead")%>%
  summarise (n())
deaths_3

denom_3 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat != "Less than 1 year" & graft_years_cat != "1 - 2 years")%>%
  summarise (n())
denom_3

p3 = deaths_3/denom_3

# Graft years = 3 - 4 years
deaths_4 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat == "3 - 4 years" & lastknownstatus == "Dead")%>%
  summarise (n())
deaths_4

denom_4 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat != "Less than 1 year" & graft_years_cat != "1 - 2 years" & graft_years_cat != "2 - 3 years")%>%
  summarise (n())
denom_4

p4 = deaths_4/denom_4

# Graft years = 4 - 5 years
deaths_5 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat == "4 - 5 years" & lastknownstatus == "Dead")%>%
  summarise (n())
deaths_5

denom_5 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat != "Less than 1 year" & graft_years_cat != "1 - 2 years" & graft_years_cat != "2 - 3 years" & graft_years_cat != "3 - 4 years")%>%
  summarise (n())
denom_5

p5 = deaths_5/denom_5

# Graft years = More than 5 years
deaths_6 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat == "More than 5 years" & lastknownstatus == "Dead")%>%
  summarise (n())
deaths_6

denom_6 = transplant_deceased_cl_selected %>%
  filter(kdpi < 25 & graftno == "1" & graft_years_cat != "Less than 1 year" & graft_years_cat != "1 - 2 years" & graft_years_cat != "2 - 3 years" & graft_years_cat != "3 - 4 years" & graft_years_cat != "4 - 5 years")%>%
  summarise (n())
denom_6

p6 = (deaths_6/denom_6)/5


p1
p2
p3
p4
p5
p6


################ Calculating graft failure rates according to number of years since transplant & different KDPI levels

### Pre-processing
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, graft_years_1 = graft_years + 0.001) # The survreg function in R does not allow time = 0


# Kaplan-Meier plot
plot_kdpi_25_G = plot(survival::survfit(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ 1, data = kdpi_25), 
                      main = "KDPI < 25",
                      xlab = "Years", 
                      ylab = "Overall graft survival probability",
                      ylim=c(0.6,1))

plot_kdpi_50_G = plot(survival::survfit(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ 1, data = kdpi_50), 
                      main = "KDPI 25-49",
                      xlab = "Years", 
                      ylab = "Overall graft survival probability",
                      ylim=c(0.6,1))

plot_kdpi_75_G = plot(survival::survfit(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ 1, data = kdpi_75), 
                      xlab = "Years", 
                      ylab = "Overall graft survival probability",
                      ylim=c(0.6,1),
                      main = "KDPI 50-74")

plot_kdpi_100_G = plot(survival::survfit(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ 1, data = kdpi_100), 
                       xlab = "Years", 
                       ylab = "Overall graft survival probability",
                       ylim=c(0.6,1),
                       main = "KDPI >74")



### Weibull Regresson model

#@@@@@@@ KDPI <25
wbmod_GF_25 = survreg(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ ageattransplant, data = kdpi_25, dist='weibull') # Gender was not significant in the model.  Thus removed
#summary(wbmod_GF_25)
ConvertWeibull(wbmod_GF_25,conf.level = 0.95)

#$$ Testing Weibull vs Parametric model
kdpi_25_WB_GF = phreg(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ ageattransplant, data = kdpi_25, dist='weibull')
kdpi_25_Cox_GF = coxreg(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ ageattransplant, data = kdpi_25)
check.dist(kdpi_25_WB_GF,kdpi_25_Cox_GF)

AIC(kdpi_25_Cox_GF) # Smaller the value, better the model
AIC(wbmod_GF_25)

BIC(kdpi_25_Cox_GF)
BIC(wbmod_GF_25)

#@@@@@@@ KDPI 25-49
wbmod_GF_50 = survreg(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ ageattransplant, data = kdpi_50, dist='weibull') # Gender was not significant in the model.  Thus removed
#summary(wbmod_GF_50)
ConvertWeibull(wbmod_GF_50,conf.level = 0.95)

#@@@@@@@ KDPI 50-74
wbmod_GF_75 = survreg(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ ageattransplant, data = kdpi_75, dist='weibull') # Gender was not significant in the model.  Thus removed
#summary(wbmod_GF_75)
ConvertWeibull(wbmod_GF_75,conf.level = 0.95)

#@@@@@@@ KDPI 50-74
wbmod_GF_100 = survreg(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ ageattransplant, data = kdpi_100, dist='weibull') # Gender was not significant in the model.  Thus removed
#summary(wbmod_GF_100)
ConvertWeibull(wbmod_GF_100,conf.level = 0.95)


################ Calculating mortality rates according to number of years since transplant after re-transplant

re_transplant = filter(transplant_deceased_cl_selected, graftno > 1)
wbmod_reT = survreg(Surv(surv_years_1, lastknownstatus == "Dead") ~ ageattransplant, data = re_transplant, dist='weibull') # Gender was not significant in the model.  Thus removed
#summary(wbmod_GF_25)
ConvertWeibull(wbmod_reT,conf.level = 0.95)



# Creating CSV files to be used in STATA (Expenential model)
write.table(kdpi_25, file = "kdpi_25.csv", sep = ",")
write.table(kdpi_50, file = "kdpi_50.csv", sep = ",")
write.table(kdpi_75, file = "kdpi_75.csv", sep = ",")
write.table(kdpi_100, file = "kdpi_100.csv", sep = ",")





