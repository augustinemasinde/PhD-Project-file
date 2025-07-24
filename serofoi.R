library(serosim)
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(reshape2)
library(readxl)
library(lubridate)
library(dplyr)
library(serofoi)
library(here)

## constant force of infection model(endemic model)
#The endemic constant model is a simple model that describe the seroprevalence of chikungunya 
#disease within a population as along term transmission. The rate of transmission is constant 
#over time and seroprevalence is a cumulative process increasing monotonically with age
chikdata <- read_excel(here("data", "chikungunya_data_Uganda.xlsx"))
#select age column
foi_df <- chikdata %>% select(Age_Yrs) %>% rename(year = Age_Yrs)
#add force of infection- how the disease spreads
foi_df$foi <- c(rep(0.02,nrow(foi_df)))
#Convert age in years to actual birth year
foi_df$year <- as.numeric(format(Sys.Date(), "%Y")) - foi_df$year
foi_df <- foi_df %>% arrange(year)
foi_df <- foi_df %>% filter(year != 1891)




#acutual survey_features
survey_features<- data.frame(
 age_min = c(1,6,12,19,50,55),
 age_max = c(5,11,18,49,54,97),
  n_sample = c(118,149,144,920,51,158)
)

#simulate a serosurvey following a stepwise decreasing FOI(lambda)
#The model assumes that the FOI is a piecewise constant over age intervals but the values
#decrease in each successive interval- immunity buildup
serosurvey_constant <- simulate_serosurvey(
  "time",
  foi = foi_df,
  survey_features
 ) %>% mutate(survey_year = 2025)

#simulated dataset contains 1541 samples of individuals between 1 and 140 years

#The code below shows how to implement the constant model to the simulated serosurvey
seromodel_constant <- fit_seromodel(
  serosurvey = serosurvey_constant1,
  model_type = "constant",
  iter = 800
)

#Code to generate the plot for the model
plot_seromodel(
  seromodel_constant,
  serosurvey = serosurvey_constant,
  foi_df = foi_df,
  size_text = 6
)
