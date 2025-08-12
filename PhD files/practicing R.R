library(ggplot2)
library(tidyverse)
pacman::p_load(
  rio,          # File import
  here,         # File locator
  skimr,        # get overview of data
  tidyverse,    # data management + ggplot2 graphics 
  gtsummary,    # summary statistics and tests
  rstatix,      # summary statistics and statistical tests
  janitor,      # adding totals and percents to tables
  scales,       # easily convert proportions to percents  
  flextable, # converting tables to pretty images
  readxl      # reading excel files
)
chikdata <- read_excel(here("Data", "chikungunya_data_Uganda.xlsx"))
skim(chikdata)



library(serofoi)
# Loading and preparing data for modelling
data(chagas2012)

# Model implementation
model_constant <- fit_seromodel(
  serosurvey = chagas2012,
  model_type = "constant"
)
# Visualisation
plot_seromodel(
  model_constant,
  serosurvey = chagas2012,
  bin_serosurvey = TRUE,
  size_text = 6
)



#Constant force of infection model
max_age <- 80
foi_constant <- data.frame(
  age = seq(1, max_age, 1),
  foi = rep(0.05, max_age)
)

ggplot(foi_constant, aes(x = age, y = foi)) +
  geom_line() + theme_bw()

#
n_sample <- 15
survey_features <- data.frame(
  age_min = seq(1, max_age, 1),
  age_max = seq(1, max_age, 1),
  n_sample = rep(n_sample, max_age))
ggplot(survey_features, aes(x = age_min, y = n_sample)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "Age (years)", y = "Sample size") +
  scale_x_continuous(breaks = seq(0, max_age, 10)) +
  scale_y_continuous(breaks = seq(0, n_sample, 5))

set.seed(1234)
serosurvey_constant <- simulate_serosurvey(
  model = "age",
  foi = foi_constant,
  survey_features = survey_features
) |>
  mutate(model_type = "constant FoI")

glimpse(serosurvey_constant)

#age varying force of infection model
ages <- seq(1, max_age, 1)
foi_age_varying <- data.frame(
  age = ages,
  foi = 0.1 * exp(-0.1 * ages)
)

ggplot(foi_age_varying, aes(x = age, y = foi)) +
  geom_line()



# Hypothetical FOI function for chikungunya
# age: age in years
# month: month of the year (1–12)
# a0: baseline FOI at age 0
# age_decay: rate of FOI decline with age
# seasonal_amp: amplitude of seasonal variation (0 to 1)
# seasonal_peak: month where FOI peaks (1 = January, 7 = July, etc.)

chik_foi <- function(age, month, 
                     a0 = 0.2,         # baseline FOI (per year) at age 0
                     age_decay = 0.1,  # decay rate per year of age
                     seasonal_amp = 0.5, # 0.5 = 50% variation from baseline
                     seasonal_peak = 5   # peak month (e.g. 8 = August)
) {
  # Age effect (exponential decay)
  age_component <- a0 * exp(-age_decay * age)
  
  # Seasonal effect using sine wave
  # Shift peak month to align with seasonal_peak
  seasonal_component <- 1 + seasonal_amp * 
    sin( 2 * pi * (month - seasonal_peak) / 12 )
  
  # Combined FOI
  foi <- age_component * seasonal_component
  return(foi)
}

# Example: FOI for a 5-year-old in August
chik_foi(age = 5, month = 5)

# Example: FOI curve for ages 0–50 in peak month
ages <- 0:50
foi_peak <- chik_foi(age = ages, month = 5)
plot(ages, foi_peak, type = "l", xlab = "Age (years)", ylab = "FOI (per year)",
     main = "Hypothetical Chikungunya FOI by Age at Peak Season")

# Example: FOI across months for a 10-year-old
months <- 1:12
foi_year <- chik_foi(age = 10, month = months)
plot(months, foi_year, type = "o", xlab = "Month", ylab = "FOI (per year)",
     main = "Hypothetical Chikungunya FOI Seasonality (Age 10)")



serosurvey_age_dep <- simulate_serosurvey(
  model = "age",
  foi = foi_age_varying,
  survey_features = survey_features
) |>
  mutate(model_type = "age-dependent FoI")
glimpse(serosurvey_age_dep)



foi_time <- c(rep(0, 40), rep(1, 40))
foi_df_time <- data.frame(
  year = seq(1946, 2025, 1),
  foi = foi_time
)

#data analysis-pivoting
#import data
count_data <- import(here("Data","malaria_facility_count_data.rds"))

# import your dataset
linelist <- import(here("Data", "linelist_cleaned.rds"))

ggplot(count_data) +
  geom_col(aes(x = data_date, y = malaria_tot), width = 1)


# Pivoting data
df_long <- 
  count_data %>% 
  pivot_longer(
    cols = starts_with("malaria_"),
    names_to = "age_group",
    values_to = "counts"
  )

df_long

df_long



data1 <- data.frame(index = c(1, 2, 3, 4),
                    name = c("A", "B", "C", "D"),
                    maths = c(10, 20, 30, 40),
                    english = c(15, 25, 35, 45),
                    science = c(20, 30, 40, 50))
data_long <- data1 %>%
  pivot_longer(
    cols = c(maths, english, science),
    names_to = "subject",
    values_to = "score"
  )
data_long


ggplot(data = data_long) +
  geom_col(
    mapping = aes(x = name, y = score, fill = subject),
    width = .5
  ) + theme_bw()
#To see the groups and the number of rows in each group, pass the grouped data to tally()
  data_group <- data_long %>% group_by(subject, name) %>% tally() 

#To see just the unique groups without counts you can pass to group_keys()
  data_group_keys <- data_long %>% group_by(name) %>% group_keys()
  
data2 <- data_long %>% group_by(pass = ifelse(score >= 20, "pass", "fail")) %>% 
  tally(sort = TRUE)
  
  
  
dall_by_outcome <- linelist %>% 
  group_by(outcome)  
  

# Grouped by outcome
by_outcome <- linelist %>% 
  group_by(outcome)

# Add grouping by gender in addition
by_outcome_gender <- by_outcome %>% 
  group_by(gender, .add = TRUE)

linelist %>% 
  group_by(outcome, gender) %>% 
  tally() %>% 
  ungroup()

linelist %>% 
  group_by(outcome, gender) %>% 
  tally()




  