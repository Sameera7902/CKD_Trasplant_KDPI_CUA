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
library(forcats)
#install.packages("data.table")
library(data.table)



#### Import the dialysis data base
setwd("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA")
dialysis = read.csv("Dialysis_data.csv")


# Select the important columns
dialysis = dialysis %>%
  select(1, 2, 14, 34, 35)


# Create death column (1 = Death; 0 = Alive)
dialysis = dialysis %>%
  mutate(death = ifelse(ageatdeath >0, "Death", "Alive"))

dialysis$death = ifelse(is.na(dialysis$death), "Alive", dialysis$death)

dialysis %>%
  summarise(prop = mean(death == "Death"))


# Varibale convertion
dialysis$treatmentdate = as.character(dialysis$treatmentdate)


# Parse as dates
dialysis$treatmentdate = parse_date(dialysis$treatmentdate, format = "%d %b %Y")

# Create curretn age column
dialysis = dialysis %>%
  mutate(duration = (as.Date("2017-12-31") - treatmentdate)/365)

dialysis = dialysis %>%
  mutate(cur_age = ageatRRTstart + duration)

dialysis$cur_age = round(dialysis$cur_age, digits = 0)


# Create age categories
dialysis = mutate(dialysis, cur_age_cat = case_when(
  cur_age <=39 ~ "Less than 40",
  cur_age >= 40 & cur_age <= 44 ~ "40 - 44",
  cur_age >= 45 & cur_age <= 49 ~ "45 - 49",
  cur_age >= 50 & cur_age <= 54 ~ "50 - 54",
  cur_age >= 55 & cur_age <= 59 ~ "55 - 59",
  cur_age >= 60 & cur_age <= 64 ~ "60 - 64",
  cur_age >= 65 & cur_age <= 69 ~ "65 - 69",
  cur_age >= 70 ~ "More than 69",
)
)


dialysis = mutate(dialysis, cur_age_cat = fct_relevel(cur_age_cat,
                                                                                                      "Less than 40",
                                                                                                      "40 - 44",
                                                                                                      "45 - 49",
                                                                                                      "50 - 54",
                                                                                                      "55 - 59",
                                                                                                      "60 - 64",
                                                                                                      "65 - 69",
                                                                                                      "More than 69",
))

# save the data
save(dialysis, file='Dialysis_processedData.RData')



