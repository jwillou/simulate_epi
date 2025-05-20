setwd("/Users/jannawilloughby/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/simulate_epi/")
#setwd("/Users/jrw0107/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/simulate_epi/")
directory = getwd()
outdir    = paste(directory, "/output/", sep = "")                 # directory to save model output  
source(paste(directory, "/source/FunctionSourcer.R", sep = ''))    # source functions and set source directory

# parameter values
t_max  = 250
env_A  = 0.5
env_B  = 1
env_sd = 0.1
n_pop  = 100
drift_C = 10

# allow shifts
start_E_A = 0.2
start_E_B = 1

# heritability and migration values
compare_h_vals = c(0.50, 0.75, 0.90)
compare_m_vals = c(0.01, 0.1, 0.2)

# run equations and plot
par(mfrow = c(length(compare_h_vals), length(compare_m_vals)), mar = c(4, 4, 2, 1))
for (h_val in compare_h_vals) {
  for (m_val in compare_m_vals) {
    simple_pred    = Theory_simple(h_val, start_E_A, env_A, t_max)
    migration_pred = Theory_migration(h_val, m_val, start_E_A, start_E_B, env_A, env_B, t_max)
    full_pred      = Theory_full(h_val, m_val, start_E_A, start_E_B, env_A, env_B, env_sd, n_pop, t_max, drift_C = 1)
    WP_pred        = Theory_WP(N = n_pop, generations = t_max, mu = 1e-5, rho = 1 - h_val, initial_freqs = c(A_M = 0.10, A_U = 0.40, a_M = 0.10, a_U = 0.40))
    accum_pred_df  = Theory_accumulative_migration(h_val, m_val, start_E_A, env_B, env_A, env_sd, n_pop, t_max, drift_C = 1)
    
    # Plot setup
    plot(-100, -100, type = "l", xlim = c(0, t_max), ylim = c(0.2, 1), xlab = "Generation", ylab = "Mean Epigenetic Value", main = paste("h =", h_val, ", m =", m_val))
    
    # Theory curves

    lines(1:t_max, migration_pred, col = "darkorange", lwd = 1, lty = 1)
    lines(1:t_max, simple_pred, col = "firebrick2", lwd = 1, lty = 1)
    lines(1:t_max, accum_pred_df$Mean_E, col = "black", lwd = 1, lty = 1)  
    lines(WP_pred$Generation, WP_pred$mean_epi, col = "navy", lwd = 1, lty = 1)
    lines(full_pred$Generation, full_pred$Mean_E, col = "forestgreen", lwd = 1, lty = 1)
  }
}

# Reset layout
par(mfrow = c(1, 1))
