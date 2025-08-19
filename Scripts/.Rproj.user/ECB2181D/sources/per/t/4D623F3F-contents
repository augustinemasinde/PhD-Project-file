# Load necessary libraries
library(readxl)
library(serofoi)
library(dplyr)
library(here)
chikdata <- read_excel(here("Data", "chikungunya_data_Uganda.xlsx"))
chikdata <- chikdata %>%  select(UniqueKey, Year, Age_Yrs,IgM_CHIK)
chikdata <- chikdata %>%
  mutate(
    IgM_CHIK = recode(IgM_CHIK, "Nengative" = "Negative", "NA" = NA_character_)
  ) %>%
  filter(!is.na(IgM_CHIK))
chikdata <- chikdata %>%
  mutate(
    age_cate = case_when(
      Age_Yrs >= 1 & Age_Yrs <= 11 ~ "1-11",
      Age_Yrs >= 12 & Age_Yrs <= 18 ~ "12-18",
      Age_Yrs >= 19 & Age_Yrs <= 49 ~ "19-49",
      Age_Yrs >= 50 & Age_Yrs <= 64 ~ "50-64",
      Age_Yrs >= 65 & Age_Yrs <= 97 ~ "65-97",
      TRUE ~ NA_character_  # catches anything outside range or NA
    )
  )

##SEROFOI PACKAGE
#Serofoi package is a tool used to estimate the force of infection(FOI) from a population based serosurvey data.
# FOI quantifies the risk of susceptible individuals becoming infected over a period of time
#We gain insights into the transmission of chikungunya dynamics in uganda, a sentinel hospital based survey data
#collected between 2028 to 2024
#Serofoi is used to compare three potential scenarios of chikungunya transmission, constant endemic, 
#epidemic slow and epidemic fast

#chikungunya transmission data in 2019
chikdata2019 <- data.frame(
  survey_year = c(2019,2019,2019,2019,2019),
  n_sample = c(115, 55,169, 23,9),
  n_seropositive = c(12,3,10,1,1),
  age_min = c(1,12,19,50,65),
  age_max = c(11,18,49,64,97)
)

chikdata2019 <- chikdata2019 %>%
  mutate(n_sample = as.integer(n_sample))

# Implementation of the models

#constant endemic model2019
seromodel_constant <- fit_seromodel(
  serosurvey = chikdata2019,
  model_type = "constant",
  iter = 1000
)

#Epidemic slow 2019
foi_index <- get_foi_index(chikdata2019, group_size = 5, model_type = "time")
seromodel_time <- fit_seromodel(
  serosurvey = chikdata2019,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  foi_index = foi_index,
  iter = 2500
)

#Epidemic fast 2019
foi_index <- get_foi_index(chikdata2019, group_size = 5, model_type = "time")
seromodel_log_time <- fit_seromodel(
  serosurvey = chikdata2019,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000
)


# Visualisation of the results
plot_constant <- plot_seromodel(
  seromodel = seromodel_constant,
  serosurvey = chikdata2019,
  foi_max = 0.07,
  size_text = 6
)
plot_time <- plot_seromodel(
  seromodel = seromodel_time,
  serosurvey = chikdata2019,
  foi_max = 0.07,
  size_text = 6
)
plot_log_time <- plot_seromodel(
  seromodel = seromodel_log_time,
  serosurvey = chikdata2019,
  foi_max = 0.07,
  size_text = 6
)

cowplot::plot_grid(plot_constant, plot_time, plot_log_time, ncol = 3)



#chikungunya transmission data  in 2020
chikdata2020 <- data.frame(
  survey_year = c(2020,2020,2020,2020,2020),
  n_sample = c(77, 31,292, 49,21),
  n_seropositive = c(17,9,57,13,2),
  age_min = c(1,12,19,50,65),
  age_max = c(11,18,49,64,97)
)

chikdata2020 <- chikdata2020 %>%
  mutate(n_sample = as.integer(n_sample))

# Implementation of the models

#constant endemic 2020
seromodel_constant <- fit_seromodel(
  serosurvey = chikdata2020,
  model_type = "constant",
  iter = 1000
)

#Epidemic slow 2020
foi_index <- get_foi_index(chikdata2020, group_size = 5, model_type = "time")
seromodel_time <- fit_seromodel(
  serosurvey = chikdata2020,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  foi_index = foi_index,
  iter = 2500
)

#Epidemic fast 2020
foi_index <- get_foi_index(chikdata2020, group_size = 5, model_type = "time")
seromodel_log_time <- fit_seromodel(
  serosurvey = chikdata2020,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000
)


# Visualisation of the results
plot_constant <- plot_seromodel(
  seromodel = seromodel_constant,
  serosurvey = chikdata2020,
  foi_max = 0.07,
  size_text = 6
)
plot_time <- plot_seromodel(
  seromodel = seromodel_time,
  serosurvey = chikdata2020,
  foi_max = 0.07,
  size_text = 6
)
plot_log_time <- plot_seromodel(
  seromodel = seromodel_log_time,
  serosurvey = chikdata2020,
  foi_max = 0.07,
  size_text = 6
)

cowplot::plot_grid(plot_constant, plot_time, plot_log_time, ncol = 3)


#chikungunya transmission data in 2021
chikdata2021 <- data.frame(
  survey_year = c(2021,2021,2021,2021,2021),
  n_sample = c(51, 42,345, 48,22),
  n_seropositive = c(3,5,24,2,2),
  age_min = c(1,12,19,50,65),
  age_max = c(11,18,49,64,97)
)

chikdata2021 <- chikdata2021 %>%
  mutate(n_sample = as.integer(n_sample))

# Implementation of the models

# Constant endemic 2021
seromodel_constant <- fit_seromodel(
  serosurvey = chikdata2021,
  model_type = "constant",
  iter = 1000
)


# Epidemic slow 2021
foi_index <- get_foi_index(chikdata2021, group_size = 5, model_type = "time")
seromodel_time <- fit_seromodel(
  serosurvey = chikdata2021,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  foi_index = foi_index,
  iter = 2500
)

#Epidemic fast 2021
foi_index <- get_foi_index(chikdata2021, group_size = 5, model_type = "time")
seromodel_log_time <- fit_seromodel(
  serosurvey = chikdata2021,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000
)


# Visualisation of the results
plot_constant <- plot_seromodel(
  seromodel = seromodel_constant,
  serosurvey = chikdata2021,
  foi_max = 0.07,
  size_text = 6
)
plot_time <- plot_seromodel(
  seromodel = seromodel_time,
  serosurvey = chikdata2021,
  foi_max = 0.07,
  size_text = 6
)
plot_log_time <- plot_seromodel(
  seromodel = seromodel_log_time,
  serosurvey = chikdata2021,
  foi_max = 0.07,
  size_text = 6
)

cowplot::plot_grid(plot_constant, plot_time, plot_log_time, ncol = 3)


#chikungunya transmission data in 2022
chikdata2022 <- data.frame(
  survey_year = c(2022,2022,2022,2022,2022),
  n_sample = c(6, 7,56, 6,6),
  n_seropositive = c(0,1,3,1,0),
  age_min = c(1,12,19,50,65),
  age_max = c(11,18,49,64,97)
)

chikdata2022 <- chikdata2022 %>%
  mutate(n_sample = as.integer(n_sample))

# Implementation of the models

#Constant endemic 2022
seromodel_constant <- fit_seromodel(
  serosurvey = chikdata2022,
  model_type = "constant",
  iter = 1000
)


#Epidemic slow 2022
foi_index <- get_foi_index(chikdata2022, group_size = 5, model_type = "time")
seromodel_time <- fit_seromodel(
  serosurvey = chikdata2022,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  foi_index = foi_index,
  iter = 2500
)

#Epidemic fast 2022
foi_index <- get_foi_index(chikdata2022, group_size = 5, model_type = "time")
seromodel_log_time <- fit_seromodel(
  serosurvey = chikdata2022,
  model_type = "time",
  foi_prior = sf_normal(0, 0.01),
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000
)


# Visualisation of the results
plot_constant <- plot_seromodel(
  seromodel = seromodel_constant,
  serosurvey = chikdata2022,
  foi_max = 0.07,
  size_text = 6
)
plot_time <- plot_seromodel(
  seromodel = seromodel_time,
  serosurvey = chikdata2022,
  foi_max = 0.07,
  size_text = 6
)
plot_log_time <- plot_seromodel(
  seromodel = seromodel_log_time,
  serosurvey = chikdata2022,
  foi_max = 0.07,
  size_text = 6
)

cowplot::plot_grid(plot_constant, plot_time, plot_log_time, ncol = 3)



