setwd("/Users/jannawilloughby/Desktop/")
library(dplyr)

# Load the saved RDS
sim_data <- readRDS("output/simulation_output_20250429_061049.rds")

# Extract everything
all_df = sim_data$all_df
agg_df = sim_data$agg_df
n_per_pop = sim_data$n_per_pop
n_gens = sim_data$n_gens
n_reps = sim_data$n_reps
h_vals = sim_data$h_vals
m_AB_vals = sim_data$m_AB
m_BA_val  = sim_data$m_BA
env_A  = sim_data$env_A
env_B  = sim_data$env_B
env_sd = sim_data$env_sd
n_loci = sim_data$n_loci