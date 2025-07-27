library(serosolver)
## Read in data
raw_dat_path <- system.file("extdata", "Fluscape_HI_data.csv", package = "serosolver")
raw_dat <- read.csv(file = raw_dat_path, stringsAsFactors = FALSE)
print(head(raw_dat))
#>   Age HI.H3N2.1968 HI.H3N2.1975 HI.H3N2.1979 HI.H3N2.1989 HI.H3N2.1995
#> 1  75           80           40           40           80          160
#> 2  35           20           80          160           40           80
#> 3  71           80           40           20           20           40
#> 4  65           80           40           40           20           40
#> 5  64          160           80           40           10           40
#> 6  33           40           20          160           80           80
#>   HI.H3N2.2002 HI.H3N2.2003 HI.H3N2.2005 HI.H3N2.2008
#> 1          160           40           80           40
#> 2           20           10            0            0
#> 3           80           20           10            0
#> 4           20            0            0            0
#> 5           40            0           20           20
#> 6          160           40           40           20

## Add indexing column for each individual
raw_dat$individual <- 1:nrow(raw_dat)

## Convert data to long format
melted_dat <- reshape2::melt(raw_dat, id.vars=c("individual","Age"),stringsAsFactors=FALSE)

## Modify column names to meet serosolver's expectations
colnames(melted_dat) <- c("individual","DOB","virus","titre")
melted_dat$virus <- as.character(melted_dat$virus)

## Extract circulation years for each virus code, which will be used 
## by serosolver as the circulation time
melted_dat$virus <- as.numeric(sapply(melted_dat$virus, function(x) strsplit(x,split = "HI.H3N2.")[[1]][2]))

## Clean and log transform the data
melted_dat <- melted_dat[complete.cases(melted_dat),]
melted_dat[melted_dat$titre == 0,"titre"] <- 5
melted_dat$titre <- log2(melted_dat$titre/5)

## Convert ages to DOB
melted_dat$DOB <- sample_year - melted_dat$DOB

## All samples taken at the same time
melted_dat$samples <- sample_year

## Add column for titre repeats, enumerating for each measurement for the same virus/sample/individual
melted_dat <- plyr::ddply(melted_dat,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x),"group"=1))

## Rename to data expected by serosolver
titre_dat <- melted_dat
print(head(titre_dat))
#>   individual  DOB virus titre samples run group
#> 1          1 1934  1968     4    2009   1     1
#> 2          1 1934  1975     3    2009   1     1
#> 3          1 1934  1979     3    2009   1     1
#> 4          1 1934  1989     4    2009   1     1
#> 5          1 1934  1995     5    2009   1     1
#> 6          1 1934  2002     5    2009   1     1


## Read in raw coordinates
antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
antigenic_coords <- read.csv(file = antigenic_coords_path, stringsAsFactors=FALSE)
print(head(antigenic_coords))
#>   Strain    X    Y
#> 1   HK68  1.8  2.4
#> 2   EN72  2.7  4.9
#> 3   VI75  7.6  6.3
#> 4   TX77  7.9  8.8
#> 5   BK79  9.6 11.0
#> 6   SI87 15.0  6.8

## Convert to form expected by serosolver
#antigenic map
antigenic_map <-example_antigenic_map


## Restrict entries to years of interest. Entries in antigenic_map determine
## the times that individual can be infected ie. the dimensions of the infection
## history matrix.
antigenic_map <- antigenic_map[antigenic_map$inf_times >= 1968 & antigenic_map$inf_times <= sample_year,]
strain_isolation_times <- unique(antigenic_map$inf_times)
print(head(strain_isolation_times))


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

## Create the posterior solving function that will be used in the MCMC framework 
model_func <- create_posterior_func(par_tab=par_tab,
                                    antibody_data = NULL,
                                    titre_dat=titre_dat,
                                    antigenic_map=antigenic_map,
                                    version=prior_version) # function in posteriors.R
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
    start_inf <- setup_infection_histories_titre(titre_dat, strain_isolation_times, 
                                                 space=3,titre_cutoff=4)
    start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
  }
  
  
  res <- serosolver(par_tab = start_tab, 
                    titre_dat = titre_dat,
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


