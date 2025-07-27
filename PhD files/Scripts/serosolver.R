install.packages("doMC")
library(serosolver)
library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(tidyr)
library(doParallel)
library(foreach)
library(ggpubr)
library(bayesplot)
library(viridis)
library(ggthemes)
library(cowplot)
library(grid)
library(gridExtra)
library(doRNG)
library(serosim)
library(tidyverse)
library(readxl)
library(doMC)
library(here)
set.seed(0)

# Setup -------------------------------------------------------------------
serosolver <- FALSE

## Filename prefix for all saved outputs
filename <- "chikungunya_data_test"
filenames <- paste0(filename, "_",1:5)
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
registerDoMC(cores=5)
stopCluster(cl)
foreach(i = 1:5) %dopar% { Sys.getpid() }

#import data
chikdata <-read_excel(here("Data", "chikungunya_data_Uganda.xlsx"))

set.seed(123) 

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
chikdata<- chikdata %>%  select(UniqueKey, Age_Yrs, Year,IgM_CHIK, titre) %>% rename(DOB = Age_Yrs)


#recode Nengative to Negative
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "Nengative"] <- "Negative"

# Recode "NA" string to actual NA
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "NA"] <- NA

# Drop all real NA values
chikdata <- chikdata[!is.na(chikdata$IgM_CHIK), ]

chikdata <- chikdata %>% dplyr::rename(individual = UniqueKey)
chikdata <- chikdata %>% select(individual,DOB,titre,Year)
chikdata <- chikdata[!is.na(chikdata$DOB), ]
chikdata$individual <- as.integer(chikdata$individual)
chikdata <- chikdata[!chikdata$Year %in% c(2024, 2025), ]


set.seed(123)
##cross sectional analysis for 2019
chikdata2019 <- chikdata[chikdata$Year == 2019, ]
chikdata2019$virus <- c(rep(2019, nrow(chikdata2019)))
chikdata2019$titre <- log2(chikdata2019$titre)
chikdata2019$samples <- 2019
chikdata2019 <- chikdata2019[!is.na(chikdata2019$titre), ]
chikdata2019$DOB <- as.integer(chikdata2019$DOB)

# plot the data
ggplot(chikdata2019, aes(x = titre)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "Distribution of Chikungunya Antibody Titres2019",
       x = "Antibody Titre (log2 scale)",
       y = "Frequency") +
  theme_minimal()

sample_year <- 2019 # Year of the serosurvey

chikdata2019$DOB <- sample_year - chikdata2019$DOB # Convert DOB to age in years at the time of the survey

chikdata2019 <- chikdata2019 %>%  select(individual,DOB,virus, titre,samples)

## Add column for titre repeats, enumerating for each measurement for the same virus/sample/individual
chikdata2019 <- plyr::ddply(chikdata2019,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x),"group"=1))
chikdata2019$individual <- 1:nrow(chikdata2019)


#antigenic map
antigenic_map <-example_antigenic_map
write.csv(antigenic_map, "antigenic_map.csv")

sample_year <- 2019 # Year of the serosurvey
#strain isolation times
antigenic_map <- antigenic_map[antigenic_map$inf_times >= 1945 & antigenic_map$inf_times <= sample_year,]
strain_isolation_times <- unique(antigenic_map$inf_times)


par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
par_tab <- read.csv(file = par_tab_path, stringsAsFactors=FALSE)

## Set parameters for Beta prior on infection histories
beta_pars <- find_beta_prior_mode(0.15,4)
par_tab[par_tab$names == "alpha","values"] <- beta_pars$alpha
par_tab[par_tab$names == "beta","values"] <- beta_pars$beta

## Maximum recordable log titre in these data is 8
par_tab[par_tab$names == "MAX_TITRE","values"] <- 8

## Remove phi parameters, as these are integrated out under prior version 2
par_tab <- par_tab[par_tab$names != "phi",]

## Fix all short term parameters to 0
par_tab[par_tab$names %in% c("mu_short","sigma2","wane"),"fixed"] <- 1 # mu_short, waning and sigma2 are fixed
par_tab[par_tab$names %in% c("mu_short","sigma2","wane"),"values"] <- 0 # set these values to 0




## Distinct filename for each chain
no_chains <- 5
filenames <- paste0(filename, "_",1:no_chains)
chain_path <- sub("par_tab_base.csv","",par_tab_path)
chain_path_real <- paste0(chain_path, "cs2_real/")
chain_path_sim <- paste0(chain_path, "cs2_sim/")


antibody_data<- antibody_data %>% rename(sample_time = samples)
## Create the posterior solving function that will be used in the MCMC framework 
model_func <- create_posterior_func(par_tab=par_tab,
                                    antibody_data = antibody_data,
                                    function_type = 1,
                                    data_type = 1,
                                    n_alive = NULL,
                                    biomarker_groups_weights = NULL,
                                    titre_dat=chikdata2019,
                                    antigenic_map= antigenic_map,
                                    version=1) # function in posteriors.R
#> Creating posterior solving function...
#> 
## Generate results in parallel
res <- foreach(x = filenames, .packages = c('serosolver','data.table','plyr')) %dopar% {
  ## Not all random starting conditions return finite likelihood, so for each chain generate random
  ## conditions until we get one with a finite likelihood
  start_prob <- -Inf
  while(!is.finite(start_prob)){
    ## Generating starting antibody kinetics parameters
    start_tab <- generate_start_tab(par_tab)
    
    ## Generate starting infection history
    start_inf <- setup_infection_histories_titre(chikdata2019, strain_isolation_times, 
                                                 space=3,titre_cutoff=4)
    start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
  }
  
  
  res <- serosolver(par_tab = start_tab, 
                    titre_dat = chikdata2019,
                    antigenic_map = antigenic_map,
                    start_inf_hist = start_inf, 
                    mcmc_pars = c("iterations"=500000,"adaptive_iterations"=100000,"thin"=1000,
                                  "thin_inf_hist"=1000,"save_block"=1000,
                                  "proposal_inf_hist_time_prop"=1, "proposal_inf_hist_indiv_prop"=1,
                                  "proposal_inf_hist_group_swap_ratio"=0.8, "proposal_inf_hist_group_swap_prop"=1),
                    filename = paste0(chain_path_real,x), 
                    CREATE_POSTERIOR_FUNC = create_posterior_func, 
                    version = prior_version)
}




