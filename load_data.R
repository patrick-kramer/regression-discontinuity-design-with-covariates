## Required Libraries
library(haven)

## Load Data
data <- read_dta("./data/Card_analysis_dataset.dta")

## Restrict to sub-sample of workers with a positive unemployment spell
data <- data[data$dunempl5>0,]
attach(data)

## Assign to those people who had no job before, a duration of the previous job of zero
ind <- which(last_job==0)
last_duration[ind] <- 0

# Basic Model
Xpaper_basis=cbind(female,married,austrian,bluecollar,age, age2,lwage,lwage2,
                   endmo_dum2,endmo_dum3,endmo_dum4, endmo_dum5,endmo_dum6,endmo_dum7,
                   endmo_dum8,endmo_dum9,endmo_dum10, endmo_dum11,endmo_dum12,
                   endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                   endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15, 
                   endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21)

## Full Model
Xpaper_extended <- cbind(Xpaper_basis,firms,experience,exper2,
                         last_job,last_bluec,last_posnonedur,
                         last_recall,last_noneduration,last_breaks,high_ed,
                         iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                         reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6)

## Outcome
Y <- wage_change

## Running Variable
X <- dten1
