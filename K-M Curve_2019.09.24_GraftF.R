# CKD - KDPI cost effectiveness analysis
# Kaplan-curves
# September 2019


# get the data
# Office Comp
setwd("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA")
load("C:/Users/n10075283/Google Drive/PhD work - Sameera/00_PhD Work/HE Work/PhD Project/TreeAge Analysis/CKD_Trasplant_KDPI_CUA/HE_Sc1_transplant.RData")


# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(survival)
library(ggplot2)


### Creating Survival years column
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, surv_years_1 = surv_years + 0.001) # The survreg function in R does not allow time = 0
transplant_deceased_cl_selected = mutate(transplant_deceased_cl_selected, graft_years_1 = graft_years + 0.001) # The survreg function in R does not allow time = 0


# Data processing for K-M curve - Graft failure
cxmod = coxph(Surv(graft_years_1, endtransplantcode == "Graft failure") ~ KDPI_cat, data = transplant_deceased_cl_selected)

newdat <- expand.grid(
  KDPI_cat = levels(transplant_deceased_cl_selected$KDPI_cat)
)

cxsf <- survfit(cxmod, data = transplant_deceased_cl_selected, newdata = newdat, conf.type = "none")
str(cxsf)


surv_cxmod0 <- surv_summary(cxsf)
head(surv_cxmod0)

surv_cxmod <- cbind(surv_cxmod0,
                    newdat[as.character(surv_cxmod0$strata), ])

surv_cxmod$strata = revalue(surv_cxmod$strata, c("1" = "KDPI < 25", "2" = "KDPI 25-49", "3" = "KDPI 50-74", "4" = "KDPI > 74"))





# K-M Curve
p1 = ggsurvplot_df(surv_cxmod, color = "strata",
              legend.title = NULL, censor = FALSE)
p1 = p1 + labs(x = "Graft survival time (Years)", y= "Graft survival probability", colour = "KDPI level")
p1 = p1 + theme(legend.position = c(0.2, 0.42))
p1 = p1 + ylim(c(0.6,1))

p1



