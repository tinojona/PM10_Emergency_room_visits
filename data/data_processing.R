###############################################################################
### DATA PROCESSING
###############################################################################

# Objective:
# create an analysis ready data set

# libraries
library(tidyverse); library(zoo); library(lubridate)

# clean the env
rm(list = ls())


# read the data
data <- read.csv("data/emerg_ap.csv")


data_clean <- data |>
  mutate(

    # "old" individuals
    all_above64y = rowSums(across(c(all_above85y, all_7584y, all_6574y)), na.rm = TRUE),

    # 2-day MA of PM10
    pm10_MA_2 = rollmean(pm10, k = 2, align = "right", na.pad = TRUE),

    # 5-day MA of PM10
    pm10_MA_5 = rollmean(pm10, k = 5, align = "right", na.pad = TRUE),

    # 10-day MA of PM10
    pm10_MA_10 = rollmean(pm10, k = 10, align = "right", na.pad = TRUE),

    # 2-day MA of no2
    no2_MA_2 = rollmean(no2, k = 2, align = "right", na.pad = TRUE),

    # 2-day MA of no2
    o3_MA_2 = rollmean(o3, k = 2, align = "right", na.pad = TRUE),

    # 2-day MA of no2
    relhum_MA_2 = rollmean(relhum, k = 2, align = "right", na.pad = TRUE),

    # heat index based on 90th percentile
    heat_90 = ifelse(tmean > quantile(tmean, .9), 0, 1),

    # extreme heat index with 99th percentile
    heat_95 = ifelse(tmean > quantile(tmean, .95), 0, 1),

    # heat index based on 90th percentile for non heat days
    heat_90_no = ifelse(tmean > quantile(tmean, .9), 1, 0),

    # extreme heat index with 99th percentile for non heat days
    heat_95_no = ifelse(tmean > quantile(tmean, .95), 1, 0),

    date = as.Date(date),

    # get dow as factor
    dow = wday(date)
  )



write_csv(data_clean, "data/emerg_ap_CLEAN.csv")
