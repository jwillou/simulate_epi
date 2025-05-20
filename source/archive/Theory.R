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


