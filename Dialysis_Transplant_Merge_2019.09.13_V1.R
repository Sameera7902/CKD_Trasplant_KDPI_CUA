# HE_Preprocessing_Dialysis dataset
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

load("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA/Dialysis_processedData.RData")
load("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA/HE_Sc1_transplant.RData")



dialysis = data.table(dialysis)
transplant = data.table(transplant_deceased_cl_selected)

# Joining data sets
dialysis_transplant_merged = merge(dialysis, transplant, by = "id", all.x = TRUE )


save(dialysis_transplant_merged, file='Dialysis_Transplant_merged.RData')




