library(tidyverse)
library(lubridate)
library(gtools)
library(data.table)
library(pROC)
library(caret)

getwd()
source("code/scripts/001_load_data.R")              #Done
source("code/scripts/002_plot_raw_data.R")          #Done
source("code/scripts/003_logistic_regression.R")    #Done