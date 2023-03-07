######################################
## Some commonly used functions for
##  - NEO BMI Metabolite MR -
##         Analysis
## by: David Hughes
## date: January 7th 2021
######################################

#########################
### LOAD LIBRARIES
#########################
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(car)
library(kableExtra)
library(ggridges)
library(cluster)
library(pcaMethods)
library(MASS)
library(sfsmisc)
library(ivreg)


#########################
### source functions
#########################
source("functions/rntransform.R")
source("functions/normW.R")
source("functions/id_outliers.R")
source("functions/lmfit.R")
# source("functions/inspection_plot.R")
source("functions/ivregfit.R")
source("functions/obsmr.R")
source("functions/obsmr_2_long.R")
# source("functions/paired_outlier_filtering.R")
# source("functions/paired_metabolite_qc.R")
# source("functions/paired_metabolite_qc_plots.R")
# source("functions/perform_paired_metabolite_qc.R")
source("functions/est_beta_varexp.R")
source("functions/est_iqr.R")
source("functions/dh_forrest_plot.R")
source("functions/hgtest.R")
source("functions/na_iqr_outliers.R")
source("functions/ztest.R")
source("functions/covariate_lms.R")
source("functions/lipoprotein_wdata.R")
source("functions/similar_estimates.R")
source("functions/bmi_metabolite_plot.R")
source("functions/lipid_conc_plot_by_bmi_cat.R")


#########################
### LOAD Data
#########################
### parameter file
pfile = "parameters/pfile.txt"
par_data = read.table(pfile, header = FALSE, sep = "=", as.is = TRUE)

## load study R data object
cat(paste0("Load the study data Rdata file\n"))
w = which(par_data[,1] == "study_Rdata_file")
f = par_data[w,2]
load(f)

## Define mydata: the working data object
cat(paste0("Define mydata working data frame\n"))
mydata = study_data$working_data

## Covariates
cat(paste0("define covariates vector\n"))
covariables = study_data$covariables

## metids
cat(paste0("define metids vector\n"))
metids = study_data$metids

## Traits
cat(paste0("define traits vector\n"))
traits = study_data$traits

######################
## Redefine the feature
## annotation
#####################
cat(paste0("Re-Define feature annotation\n"))
## metabolite ids
fanno = study_data$feature_anno

## Match IDs
m = match(traits, fanno$new_id)
fanno = fanno[m, ]

## Redefine metabolite
fanno$metabolite = rep(metids, 3)

## Redefine ids
l = length(table( fanno$metabolite ) )
fanno$ids  = paste0(fanno$metabolite, "_", c( rep("f",l), rep("p",l) , rep("r",l)  ) )

#### update dietary state
w = which( is.na( fanno$dietary_state ) ); fanno$dietary_state[w] = "response"

##### add missing annotation
w = which(is.na(fanno$raw.label))
m = match( fanno$metabolite[w], fanno$metabolite )
fanno[w,5:10] = fanno[m,5:10]

## Redefine rownames
rownames(fanno) = fanno$ids
fanno$new_id = fanno$ids

## redefine feature annotation file
study_data$feature_anno = fanno
rm(fanno)

#############################
## Define BMI categories
#############################
cat(paste0("Define a BMI category factor\n"))
mydata$bmi_cat = sapply(mydata$bmi, function(x){
  if( x < 25){
    out = "healthy"
  } else {
    if( x >= 25 & x < 30){
      out = "overweight"
    } else {
      if(x >= 30 & x < 40){
        out = "obese"
      } else {
        if(x > 40){
          out = "severely obese"
        }
      }
    }
  }
  return(out)
})
mydata$bmi_cat = factor(mydata$bmi_cat, levels = c("healthy","overweight","obese","severely obese"))
