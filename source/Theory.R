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

# Theory for focal population under one-way migration from a static source.
# E_t+1 = inheritance from previous mean + environmental resetting toward a fixed target.
theory_simple <- function(h_epi, start_E, env_target, t_max) {
  E <- numeric(t_max)
  E[1] <- start_E
  for (t in 2:t_max) {
    E[t] <- h_epi * E[t-1] + (1 - h_epi) * env_target
  }
  return(E)
}

theory_migration <- function(h_epi, m, start_E_local, start_E_migrant, env_local, env_migrant, t_max) {
  E_local <- numeric(t_max)
  E_migrant <- numeric(t_max)
  
  E_local[1] <- start_E_local
  E_migrant[1] <- start_E_migrant
  
  for (t in 2:t_max) {
    mixed_E_local   <- (1 - m) * E_local[t-1] + m * E_migrant[t-1]
    mixed_E_migrant <- (1 - m) * E_migrant[t-1] + m * E_local[t-1]
    
    E_local[t]   <- h_epi * mixed_E_local + (1 - h_epi) * env_local
    E_migrant[t] <- h_epi * mixed_E_migrant + (1 - h_epi) * env_migrant
  }
  
  return(E_local)
}

theory_full <- function(h_epi, m, start_E_local, start_E_migrant, env_local, env_migrant, env_sd, n_pop, t_max, drift_C = 1) {
  
  mean_local <- numeric(t_max)
  mean_migrant <- numeric(t_max)
  var_local <- numeric(t_max)
  
  mean_local[1] <- start_E_local
  mean_migrant[1] <- start_E_migrant
  var_local[1] <- 0  # Assume no initial variance
  
  for (t in 2:t_max) {
    
    # Step 1: Migration
    E_after_migration_local   <- (1 - m) * mean_local[t-1] + m * mean_migrant[t-1]
    E_after_migration_migrant <- (1 - m) * mean_migrant[t-1] + m * mean_local[t-1]
    
    # Step 2: Deterministic mean
    deterministic_local <- h_epi * E_after_migration_local + (1 - h_epi) * env_local
    deterministic_migrant <- h_epi * E_after_migration_migrant + (1 - h_epi) * env_migrant
    
    # Step 3: Drift standard deviations
    drift_sd_local   <- sqrt(drift_C * mean_local[t-1] * (1 - mean_local[t-1]) / (2 * n_pop))
    drift_sd_migrant <- sqrt(drift_C * var_local[t-1] / (2 * n_pop))  # assuming same variance
    
    # Step 4: Update means with drift
    mean_local[t] <- deterministic_local + rnorm(1, mean = 0, sd = drift_sd_local)
    mean_migrant[t] <- deterministic_migrant + rnorm(1, mean = 0, sd = drift_sd_migrant)
    
    # Step 5: Update variance (local only)
    var_local[t] <- (h_epi^2) * var_local[t-1] + (1 - h_epi)^2 * env_sd^2
  }
  
  return(data.frame(
    Generation = 1:t_max,
    Mean_E = mean_local,
    Var_E = var_local,
    SD_E = sqrt(var_local)
  ))
}


# Set parameters to explore
compare_h_vals <- c(0.01, 0.25, 0.75)
compare_m_vals <- c(0.01, 0.2, 0.5)
t_max <- max(agg_df$Generation)
env_A_val <- env_A
env_B_val <- env_B
env_sd_fixed <- env_sd
n_pop <- n_per_pop
drift_C <- 1  # drift parameter for theory_full

# Storage
full_fit_stats <- data.frame()

# Plot setup
par(mfrow = c(length(compare_h_vals), length(compare_m_vals)), mar = c(4, 4, 2, 1))

for (h_val in compare_h_vals) {
  for (m_val in compare_m_vals) {
    
    sim_df <- agg_df %>%
      filter(h_epi == h_val, m == m_val) %>%
      arrange(Generation)
    
    if (nrow(sim_df) == 0) next
    
    # Starting points
    start_E_local <- mean(sim_df$epi_A.mean[sim_df$Generation == 1])
    start_E_migrant <- mean(sim_df$epi_B.mean[sim_df$Generation == 1])
    
    # Theoretical predictions
    simple_pred <- theory_simple(h_val, start_E_local, env_A_val, t_max)
    migration_pred <- theory_migration(h_val, m_val, start_E_local, start_E_migrant, env_A_val, env_B_val, t_max)
    full_pred_df <- theory_full(h_val, m_val, start_E_local, start_E_migrant, env_A_val, env_B_val, env_sd_fixed, n_pop, t_max, drift_C)
    full_pred <- full_pred_df$Mean_E  # Use mean from the full model output
    
    # RMSE for each model
    RMSE_simple <- sqrt(mean((sim_df$epi_A.mean - simple_pred)^2, na.rm = TRUE))
    RMSE_middle <- sqrt(mean((sim_df$epi_A.mean - migration_pred)^2, na.rm = TRUE))
    RMSE_full   <- sqrt(mean((sim_df$epi_A.mean - full_pred)^2, na.rm = TRUE))
    
    # Best model
    best_model <- c("Simple", "Middle", "Full")[which.min(c(RMSE_simple, RMSE_middle, RMSE_full))]
    
    full_fit_stats <- rbind(full_fit_stats,
                            data.frame(h_epi = h_val, m = m_val,
                                       RMSE_simple = RMSE_simple,
                                       RMSE_middle = RMSE_middle,
                                       RMSE_full = RMSE_full,
                                       Best_Model = best_model))
    
    # --- Plot ---
    plot(sim_df$Generation, sim_df$epi_A.mean, type = "n",
         ylim = c(0, 1), xlab = "Generation", ylab = "Epigenetic Value",
         main = paste0("h = ", h_val, ", m = ", m_val))
    
    # Simulated means and CI shading
    polygon(c(sim_df$Generation, rev(sim_df$Generation)),
            c(sim_df$epi_A.upr, rev(sim_df$epi_A.lwr)),
            col = adjustcolor("blue", alpha.f = 0.2), border = NA)
    lines(sim_df$Generation, sim_df$epi_A.mean, col = "blue", lwd = 2)
    
    # Theory lines
    lines(1:t_max, simple_pred, col = "darkorange", lwd = 2, lty = 3)  # Simple (inheritance + environment only)
    lines(1:t_max, migration_pred, col = "forestgreen", lwd = 2, lty = 3)  # Migration added
    lines(1:t_max, full_pred, col = "red", lwd = 2, lty = 3)  # Full model with drift
  }
}

# Reset plotting window
par(mfrow = c(1,1))

# Print fit results
print(full_fit_stats)

