library(tidyverse)
library(lubridate)
library(gtools)
library(data.table)
library(pROC)
library(caret)
library(R2jags)
library(rjags)

getwd()
source("code/scripts/001_load_data.R")                 #Done
source("code/scripts/002_plot_raw_data.R")             #Done
source("code/scripts/003_logistic_regression.R")       #Done
source("code/scripts/004_grid_search_occupancy.R")     #Done
source("code/scripts/005_grid_search_MS_occupancy.R")  #Done
source("code/scripts/006_HMM.R")                       #Done
source("code/scripts/007_conditional_HMM.R")           #Done