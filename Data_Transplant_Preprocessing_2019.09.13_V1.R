# HE_Preprocessing
# process the data ready for analysis - Markov model Scenario 01 (Transplant data set)
# July 2019

# Load libraries
#install.packages("haven")
library(haven)
library(dplyr)
library(readr)
#install.packages("lubridate")
library(lubridate)
#install.packages("tidyr")
library(tidyr)
library(forcats)


#### Import the deseased donor data base
setwd("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA")
transplant_deceased = read.delim("Kidney Tx_Deceased_25.03.2019.csv", sep = ",")


### Select the important columns
transplant_deceased_cl_selected = transplant_deceased %>%
  select(1:3, 17, 19, 20:27, 83:84, 94, 106, 107, 131:133)

### Check the data set
# Check the class of bmi
class (transplant_deceased_cl_selected)
# Check the dimensions of bmi
dim(transplant_deceased_cl_selected)

glimpse(transplant_deceased_cl_selected)
summary(transplant_deceased_cl_selected)

### Varibale convertion
transplant_deceased_cl_selected$endtransplantcode = as.factor(transplant_deceased_cl_selected$endtransplantcode)
transplant_deceased_cl_selected$lastknownstatus = as.factor(transplant_deceased_cl_selected$lastknownstatus)
transplant_deceased_cl_selected$gendercode = as.factor(transplant_deceased_cl_selected$gendercode)

### Parse as dates
ymd("2013-04-28")
transplant_deceased_cl_selected$transplantdate = dmy(transplant_deceased_cl_selected$transplantdate)
transplant_deceased_cl_selected$graftfailuredate = dmy(transplant_deceased_cl_selected$graftfailuredate)
transplant_deceased_cl_selected$endtransplantdate = dmy(transplant_deceased_cl_selected$endtransplantdate)
transplant_deceased_cl_selected$lastfollowupdate = dmy(transplant_deceased_cl_selected$lastfollowupdate)
transplant_deceased_cl_selected$deathdate = dmy(transplant_deceased_cl_selected$deathdate)


### Missing values
# 119 KDPI records missing

### Recording variables
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, endtransplantcode = fct_recode(endtransplantcode,
                                          "Lost to follow-up" = "K",
                                          "Graft failure" = "P",
                                          "Functioning graft" = "S",
                                          "Dead" = "Z"))

transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, lastknownstatus = fct_recode(lastknownstatus,
                                                                                                         "Lost to follow-up" = "K",
                                                                                                         "Renal recovery" = "J",
                                                                                                         "Alive" = "S",
                                                                                                         "Dead" = "Z"))


# Create age category varibale
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, final_age = (lastfollowupdate - transplantdate)/365) # Final age at follow-up is not given in the data set
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, final_age = (ageattransplant + final_age))
transplant_deceased_cl_selected$final_age = round(transplant_deceased_cl_selected$final_age, digits = 0)
transplant_deceased_cl_selected$final_age = as.numeric(transplant_deceased_cl_selected$final_age)


transplant_deceased_cl_selected %>% # To test the accuracy of the "Final_age" column
  filter(lastknownstatus == "Dead")%>%
  select(transplantdate, lastfollowupdate, ageatdeath,ageattransplant, final_age) # "Final_age" column is correct!


transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, Final_age_cat = case_when(
  final_age <=39 ~ "Less than 40",
  final_age >= 40 & final_age <= 44 ~ "40 - 44",
  final_age >= 45 & final_age <= 49 ~ "45 - 49",
  final_age >= 50 & final_age <= 54 ~ "50 - 54",
  final_age >= 55 & final_age <= 59 ~ "55 - 59",
  final_age >= 60 & final_age <= 64 ~ "60 - 64",
  final_age >= 65 & final_age <= 69 ~ "65 - 69",
  final_age >= 70 ~ "More than 69",
)
)

transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, Final_age_cat = fct_relevel(Final_age_cat,
                                                                                                      "Less than 40",
                                                                                                      "40 - 44",
                                                                                                      "45 - 49",
                                                                                                      "50 - 54",
                                                                                                      "55 - 59",
                                                                                                      "60 - 64",
                                                                                                      "65 - 69",
                                                                                                      "More than 69",
))


transplant_deceased_cl_selected$Final_age_cat = as.factor(transplant_deceased_cl_selected$Final_age_cat)


## Calculating graft failure duration category
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, graft_years = (transplantperiod/365)) # Creating "Graft survival years" variable
transplant_deceased_cl_selected$graft_years = round(transplant_deceased_cl_selected$graft_years, digits = 1)

transplant_deceased_cl_selected %>% # To test the accuracy
  filter(endtransplantcode == "Graft failure")%>%
  select(transplantdate, lastfollowupdate, ageattransplant, ageatgraftfailure, final_age, graft_years) # "graft_years" column is correct


transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, graft_years_cat = case_when(
  graft_years < 1 ~ "Less than 1 year",
  graft_years >= 1 & graft_years < 2 ~ "1 - 2 years",
  graft_years >= 2 & graft_years < 3 ~ "2 - 3 years",
  graft_years >= 3 & graft_years < 4 ~ "3 - 4 years",
  graft_years >= 4 & graft_years < 5 ~ "4 - 5 years",
  graft_years >= 5 ~ "More than 5 years",
)
)

transplant_deceased_cl_selected$graft_years_cat  = as.factor(transplant_deceased_cl_selected$graft_years_cat )

transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, graft_years_cat = fct_relevel(graft_years_cat,
                                                                                                        "Less than 1 year",
                                                                                                        "1 - 2 years",
                                                                                                        "2 - 3 years",
                                                                                                        "3 - 4 years",
                                                                                                        "4 - 5 years",
                                                                                                        "More than 5 years",
))


## Calculate patient survival years
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, surv_years = (aliveperiod/365)) # Creating "Patient survival years" variable
transplant_deceased_cl_selected$surv_years = round(transplant_deceased_cl_selected$surv_years, digits = 2)

## Creating KDPI categories
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, KDPI_cat = case_when(
  kdpi < 25 ~ "KDPI < 25",
  kdpi >= 25 & kdpi < 50 ~ "KDPI 25-49",
  kdpi >= 50 & kdpi < 75 ~ "KDPI 50-74",
  kdpi >= 75  ~ "KDPI > 74",
)
)

transplant_deceased_cl_selected$KDPI_cat  = as.factor(transplant_deceased_cl_selected$KDPI_cat )

transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, KDPI_cat = fct_relevel(KDPI_cat,
                                                                                                        "KDPI < 25",
                                                                                                        "KDPI 25-49",
                                                                                                        "KDPI 50-74",
                                                                                                        "KDPI > 74",
                                                                                                    ))


# save the data
save(transplant_deceased_cl_selected, file='HE_Sc1_transplant.RData')



