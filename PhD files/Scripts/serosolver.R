devtools::install_github("seroanalytics/serosolver", force = TRUE)
## Required for this analysis
library(serosolver)
library(reshape2)
library(foreach)
library(doParallel)
library(bayesplot)
library(coda)
library(ggplot2)
library(viridis)
library(ggpubr)
library(lubridate)
library(here)
library(readxl)
library(dplyr)
library(tidyverse)
# set up cluster
set.seed(0)
cl <- makeCluster(5)

## and the setup for parallelisation is different on a Linux or Mac:

if(Sys.info()[["sysname"]]=="Darwin" | Sys.info()[["sysname"]]=="Linux"){
  library(doMC)
  library(doRNG)
  registerDoMC(cores=5)
}else{
  registerDoParallel(cl)
}

# Set the prior version
prior_version <- 2

#import data
chikdata <- read_excel(here("Data", "chikungunya_data_Uganda.xlsx"))


chikdata$titre <- NA_real_

# Assign low titres to seronegative individuals
chikdata$titre[chikdata$IgM_CHIK == "Negative"] <- runif(
  sum(chikdata$IgM_CHIK == "Negative"), min = 0, max = 10
)

# Assigning higher titres to seropositive individuals
chikdata$titre[chikdata$IgM_CHIK == "Positive"] <- rlnorm(
  sum(chikdata$IgM_CHIK == "Positive"), meanlog = log(40), sdlog = 0.4
)

#select relevant columns
chikdata<- chikdata %>%  select(UniqueKey, Age_Yrs, Year,IgM_CHIK, titre) %>% dplyr::rename(DOB = Age_Yrs)



#recode Nengative to Negative
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "Nengative"] <- "Negative"

# Recode "NA" string to actual NA
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "NA"] <- NA

# Drop all real NA values
chikdata <- chikdata[!is.na(chikdata$IgM_CHIK), ]

chikdata <- chikdata %>% dplyr::rename(individual = UniqueKey)
chikdata <- chikdata %>% select(individual,DOB,titre)
chikdata <- chikdata[!is.na(chikdata$DOB), ]
chikdata$individual <- as.integer(chikdata$individual)
chikdata$titre <- round(chikdata$titre, 2)

filename <- "chik_2019"
resolution <- 4 # Resolution of the serosurvey, on quartely basis
sample_year <- 2019: 2022 # Year of the serosurvey

#Data preparation(create a hypothetical serosurvey data that serosolver can use)

# Create a sequence of sample years
years <- 2019:2022

# Create a sequence of sample years for each individual
chikdata <- chikdata %>%
  slice(rep(1:n(), each = 4)) %>%  
  mutate(year = rep(years, times = nrow(chikdata)))

set.seed(123)
# Randomly assign titres to each quarter
chikdata <- chikdata %>%
  group_by(individual) %>%
  mutate(
    titre = round(titre * runif(n(), 0.9, 1.1), 2)
  ) %>%
  ungroup()
# Create a virus column and assign a year of isolation
chikdata$virus <- 8077

set.seed(123)  # For reproducibility

chikdata <- chikdata %>%
  mutate(
    quarter = sample(1:4, size = n(), replace = TRUE), 
    samples = (year + (quarter - 1) / 4) * 4 + 1,   
    samples = as.integer(samples)              
  )

# Convert DOB to Date format and compute age
chikdata <- chikdata %>%
  mutate(
    DOB = as.Date(DOB),                        # convert to Date format
    birth_year = year(DOB),                    # extract year
    age = 2019 - birth_year,                   # compute age in 2019
    DOB = (birth_year * 4) + 1   # convert to sample time format
  )

indivs <- unique(chikdata$individual) #all individuals


# Subset data for indivs
chikdata <- chikdata[chikdata$individual %in% indivs,
                       c("individual","virus","titre","samples","DOB")]
chikdata$individual <- match(chikdata$individual, indivs)


# Ensure that the data is unique
chikdata <- unique(chikdata)
chikdata <- plyr::ddply(chikdata, ~individual + virus + samples,
                        function(x) cbind(x, run = 1:nrow(x), group = 1))

# Convert virus, titre, and DOB to integer
chikdata <- chikdata %>%
  mutate(
    virus = as.integer(virus),
    titre = as.integer(titre),
    DOB = as.integer(DOB)  
  )
#convert titre to log2 scale
chikdata$titre <- log2(chikdata$titre)

strain_isolation_times <- seq(sample_years[1]*resolution+1, sample_years[4]*resolution, by=1)


par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
par_tab <- read.csv(par_tab_path, stringsAsFactors=FALSE)

## Set parameters for beta and alpha to 1
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(1/3,1/3)
## Maximum recordable log titre in these data is 9
par_tab[par_tab$names == "MAX_TITRE","values"] <- 9

## Remove phi parameters, as these are integrated out under prior version 2
par_tab <- par_tab[par_tab$names != "phi",]

## Fix cross reactivity and antigenic seniority
par_tab[par_tab$names %in% c("tau","sigma1","sigma2"),"fixed"] <- 1 
## mu, tau, sigma1, and sigma2 are fixed
par_tab[par_tab$names %in% c("tau","sigma1","sigma2"),"values"] <- 0 
## set these values to 0
#possible exposure times
possible_exposure_times <- seq(1928,2022,by=1)
## Distinct filename for each chain
no_chains <- 5
filenames <- paste0(filename, "_",1:no_chains)
chain_path <- sub("par_tab_base.csv","",par_tab_path)
chain_path_real <- paste0(chain_path, "cs1_real/")
chain_path_sim <- paste0(chain_path, "cs1_sim/")

## Create the posterior solving function that will be used in the MCMC framework 
par_tab[par_tab$names == "mu_short","lower_bound"] <- 1
model_func <- create_posterior_func(par_tab=par_tab,
                                    antibody_data = NULL, # no antibody data
                                    chikdata,
                                    strain_isolation_times = strain_isolation_times,
                                    possible_exposure_times=possible_exposure_times,
                                    version=prior_version) # function in posteriors.R
#> Creating posterior solving function...
## Generate results in parallel
res <- foreach(x = filenames, .packages = c('serosolver','data.table','plyr')) %dopar% {
  ## Not all random starting conditions return finite likelihood, so for each chain generate random
  ## conditions until we get one with a finite likelihood
  start_prob <- -Inf
  while(!is.finite(start_prob)){
    ## Generating starting antibody kinetics parameters
    start_tab <- generate_start_tab(par_tab)
    
    ## Generate starting infection history
    start_inf <- setup_infection_histories_titre(chikdata, strain_isolation_times, 
                                                 space=3,titre_cutoff=4)
    start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
  }
  
  res <- serosolver(par_tab = start_tab, 
                    titre_dat = chikdata,
                    antigenic_map = NULL,
                    strain_isolation_times = strain_isolation_times,
                    start_inf_hist = start_inf, 
                    mcmc_pars = c("iterations"=2000000,"target_acceptance_rate_theta"=0.44,"target_acceptance_rate_inf_hist"=0.44,
                                  "adaptive_frequency"=1000,"thin"=1,"adaptive_iterations"=500000, 
                                  "save_block"=1000, "thin_inf_hist"=100,"proposal_inf_hist_indiv_prop"=1,
                                  "proposal_ratio"=2, "burnin"=0, "proposal_inf_hist_time_prop"=0.5, 
                                  "proposal_inf_hist_distance"=3,"proposal_inf_hist_adaptive"=1,"proposal_inf_hist_indiv_swap_ratio"=0.5,
                                  "proposal_inf_hist_group_swap_ratio"=0.5,"proposal_inf_hist_group_swap_prop"=1),
                    filename = paste0(chain_path_real,x), 
                    CREATE_POSTERIOR_FUNC = create_posterior_func, 
                    version = prior_version)
}

