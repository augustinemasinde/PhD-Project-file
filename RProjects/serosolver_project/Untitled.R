# Set seed for reproducibility
set.seed(1)

# Load required libraries
library(tidyverse)
library(devtools)

# Load serosolver package from local source
devtools::load_all("serosolver")

# Load antigenic map data
antigenic_map <- read.csv("serosolver/inst/extdata/antigenic_maps/antigenicMap_vietnam.csv")
possible_exposure_times <- antigenic_map$inf_times

# Define sampled antigens (every 2 years)
sampled_antigens <- seq(min(possible_exposure_times), max(possible_exposure_times), by = 2)

# Define sampling times and number of samples
sampling_times <- 2010:2015
n_samps <- 5

# Create attack_rates data.frame with required columns
attack_rates <- data.frame(
  population_group = 1,
  time = possible_exposure_times,
  attack_rate = runif(length(possible_exposure_times), 0.05, 0.15)
)
# Add prob_infection column if needed by package internals
attack_rates$prob_infection <- attack_rates$attack_rate

# Load parameter table
par_tab <- read.csv("serosolver/inst/extdata/par_tab_base.csv")

# Check if 'wane_maternal' parameter exists, if not add it with default values
if (!("wane_maternal" %in% par_tab$names)) {
  new_param <- data.frame(
    names = "wane_maternal",
    values = 0,
    fixed = 1,
    lower_bound = 0,
    upper_bound = 10,
    lower_start = 0,
    upper_start = 0,
    par_type = 1,
    stratification = NA,
    biomarker_group = 1
  )
  par_tab <- bind_rows(par_tab, new_param)
}

# Run the simulation
all_simulated_data <- simulate_data(
  par_tab = par_tab,
  group = 1,
  n_indiv = 50,
  possible_exposure_times = possible_exposure_times,
  measured_biomarker_ids = sampled_antigens,
  sampling_times = sampling_times,
  nsamps = n_samps,
  antigenic_map = antigenic_map,
  age_min = 10,
  age_max = 75,
  attack_rates = attack_rates,
  repeats = 2,
  data_type = c(1)
)

# Extract simulated data
antibody_data <- all_simulated_data$antibody_data
true_inf_hist <- all_simulated_data$infection_histories

# Plot antibody data
plot_antibody_data(antibody_data, possible_exposure_times, 1:4, infection_histories = true_inf_hist)

# Create data directory if doesn't exist
dir.create("serosolver/data", recursive = TRUE, showWarnings = FALSE)

# Save example data for later use
save(par_tab, file = "serosolver/data/example_par_tab.RData")
save(antibody_data, file = "serosolver/data/example_antibody_data.RData")
save(true_inf_hist, file = "serosolver/data/example_inf_hist.RData")
