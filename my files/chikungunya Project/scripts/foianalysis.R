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

##import the data
chikdata <- read_excel(here("Data", "chikungunya_data_Uganda.xlsx"))
#select age column
foi_df <- chikdata %>% select(Age_Yrs) %>% rename(year = Age_Yrs)
#add force of infection- how the disease spreads
foi_df$foi <- c(rep(0.02,nrow(foi_df)))
#omit NA
foi_df <- na.omit(foi_df)
#Convert age in years to actual birth year
foi_df$year <- as.numeric(format(Sys.Date(), "%Y")) - foi_df$year
foi_df <- foi_df %>% arrange(year)
foi_df <- foi_df %>% filter(year != 1891)




#survey features
#age categor interms age_min,age_max and n_sample 
survey_features <- chikdata %>% select(Age_Yrs)
breaks <- c(1, 5, 11, 18, 49, 54, 98)
labels <- c("1-5", "6-11", "12-18", "19-49", "50-54", "55-100")

survey_features$age_cat <- cut(
  survey_features$Age_Yrs,
  breaks = breaks,
  labels = labels,
  right = TRUE,        
  include.lowest = TRUE  
)


survey_features <- survey_features %>%
  group_by(age_cat) %>%
  mutate(n_sample = n()) %>%
  ungroup() %>% distinct()

survey_features <- survey_features %>%
  distinct(age_cat, .keep_all = TRUE)


survey_features <- survey_features %>%
  filter(!is.na(age_cat))

#acutual survey_features
survey_features1<- data.frame(
  age_min = c(1,6,12,19,50,55),
  age_max = c(5,11,18,49,54,97),
  n_sample = c(118,149,144,920,51,158)
)


#simulate a serosurvey following a stepwise decreasing FOI(lambda)
#The model assumes that the FOI is a piecewise constant over age intervals but the values
#decrease in each successive interval- immunity buildup
serosurvey_constant1 <- simulate_serosurvey(
  "time",
  foi_df,
  survey_features1
) %>% mutate(survey_year = 2025)

#simulated dataset contains 1541 samples of individuals between 1 and 140 years

#The code below shows how to implement the constant model to the simulated serosurvey
seromodel_constant <- fit_seromodel(
  serosurvey = serosurvey_constant,
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

# Time- varying models FOI models
# Time-varying model uses a forward random walk algorithm of the first chronological FOI value
#in the time-span of the of the serological survey

