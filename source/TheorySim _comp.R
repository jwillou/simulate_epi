setwd("/Users/jannawilloughby/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/simulate_epi/")
#setwd("/Users/jrw0107/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/simulate_epi/")
directory = getwd()
outdir    = paste(directory,"/output/",sep="")                    #directory to save model output  
source(paste(directory, "/source/FunctionSourcer.R", sep = ''))   #source functions and set source directory

# Load the saved RDS
sim_data = readRDS("../output/simulation_output_20250519_191149.rds")

# Extract everything
all_df    = sim_data$all_df
agg_df    = sim_data$agg_df
n_per_pop = sim_data$n_per_pop
n_gens    = sim_data$n_gens
n_reps    = sim_data$n_reps
h_vals    = sim_data$h_vals
m_AB_vals = sim_data$m_AB_vals 
m_BA_val  = sim_data$m_BA_val
env_A     = sim_data$env_A
env_B     = sim_data$env_B
env_sd    = sim_data$env_sd
n_loci    = sim_data$n_loci

# Initialize empty matrices for each model
RMSE_simple = RMSE_migration = RMSE_full = RMSE_WP = RMSE_accumulation = matrix(NA_real_, nrow = length(h_vals), ncol = length(m_AB_vals), dimnames = list(h_vals, m_AB_vals))

# RMSE helper
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
par(mfrow = c(1, 1))
# Loop over all h and m combinations
for (i in seq_along(h_vals)) {
  h = h_vals[i]
  for (j in seq_along(m_AB_vals)) {
    m = m_AB_vals[j]
    
    # Simulation subset for this combo
    sim_subset = agg_df[agg_df$h_epi == h & agg_df$m_AB == m, ]
    sim_subset = sim_subset[order(sim_subset$Generation), ]
    if (nrow(sim_subset) == 0) next
    
    sim_vec = sim_subset$epi_A.mean
    
    # Starting points for theory
    start_E_A = sim_subset$epi_A.mean[1]
    start_E_B = sim_subset$epi_B.mean[1]
    
    # Predict from theory
    pred_simple    = Theory_simple(h, start_E_A, env_A, n_gens)
    pred_migration = Theory_migration(h, m, start_E_A, start_E_B, env_A, env_B, n_gens)
    pred_full_df   = Theory_full(h, m, start_E_A, start_E_B, env_A, env_B, env_sd, n_per_pop, n_gens, drift_C = 1)
    pred_full      = pred_full_df$Mean_E
    pred_WP        = Theory_WP(N = n_per_pop, generations = n_gens, mu = 1e-5, rho = 1 - h, initial_freqs = c(A_M = 0.10, A_U = 0.40, a_M = 0.10, a_U = 0.40))$mean_epi
    pred_accum_df  = Theory_accumulative_migration(h_epi = h, m = m, start_E_A = start_E_A, env_B = env_B, env_A = env_A, env_sd = env_sd, n_pop = n_per_pop, t_max = n_gens, drift_C = 10)
    pred_accum     = pred_accum_df$Mean_E
    
    # Compute RMSE for each theory
    RMSE_simple[i, j]       = rmse(sim_vec, pred_simple)
    RMSE_migration[i, j]    = rmse(sim_vec, pred_migration)
    RMSE_full[i, j]         = rmse(sim_vec, pred_full)
    RMSE_WP[i, j]           = rmse(sim_vec, pred_WP)
    RMSE_accumulation[i, j] = rmse(sim_vec, pred_accum)
  }
}

# Display results
cat("RMSE - Theory_simple:\n");       print(round(RMSE_simple, 4))
cat("RMSE - Theory_migration:\n");    print(round(RMSE_migration, 4))
cat("RMSE - Theory_full:\n");         print(round(RMSE_full, 4))
cat("RMSE - Theory_WP:\n");           print(round(RMSE_WP, 4))
cat("RMSE - Theory_accumulation:\n"); print(round(RMSE_accumulation, 4))

filled.contour(
  x = as.numeric(colnames(RMSE_accumulation)),
  y = as.numeric(rownames(RMSE_accumulation)),
  z = t(RMSE_accumulation),
  color.palette = colorRampPalette(c("white", "orange", "red")),
  xlab = "Migration rate (m)",
  ylab = "Heritability (h_epi)",
  main = "RMSE – Theory_accumulation"
)

# compare all
rmse_list <- list(
  Simple = as.vector(RMSE_simple),
  Migration = as.vector(RMSE_migration),
  Full = as.vector(RMSE_full),
  Accum = as.vector(RMSE_accumulation),
  WP = as.vector(RMSE_WP)
)

# Compute mean and standard error for each
mean_RMSEs <- sapply(rmse_list, function(x) mean(x, na.rm = TRUE))
se_RMSEs   <- sapply(rmse_list, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Barplot
bp <- barplot(mean_RMSEs, names.arg = names(mean_RMSEs), col = c("firebrick2", "darkorange", "forestgreen", "black", "navy"), ylim = c(0, max(mean_RMSEs + se_RMSEs, na.rm = TRUE) * 1.1), ylab = "Mean RMSE",  main = "Average RMSE across scenarios")
arrows(x0 = bp, y0 = mean_RMSEs - se_RMSEs, x1 = bp, y1 = mean_RMSEs + se_RMSEs, angle = 90, code = 3, length = 0.05)


# Plot simulation + theory side-by-side 
# Choose combos to show (or leave as full vectors)
h_show  = 0.75#c(0.50, 0.75, 0.90) #h_vals
cm_show = 0.1#c(0.01, 0.1, 0.2) #m_AB_vals

# Plot simulation and theory 
par(mfrow = c(length(h_show), length(cm_show)), mar = c(4, 4, 2, 1))

for (h in h_show) {
  for (m in cm_show) {
    # Subset simulation data
    sim_df = subset(agg_df, h_epi == h & m_AB == m)
    sim_df = sim_df[order(sim_df$Generation), ]
    if (nrow(sim_df) == 0L) next
    
    #  Starting trait values for theory
    start_E_A = sim_df$epi_A.mean[1]
    start_E_B = sim_df$epi_B.mean[1]
    
    # Theory predictions
    simple_pred    = Theory_simple(h, start_E_A, env_A, n_gens)
    migration_pred = Theory_migration(h, m, start_E_A, start_E_B, env_A, env_B, n_gens)
    full_pred      = Theory_full(h, m, start_E_A, start_E_B, env_A, env_B, env_sd, n_per_pop, n_gens, drift_C = 10)
    WP_pred        = Theory_WP(N = n_per_pop, generations = n_gens, mu = 1e-5, rho = 1 - h, initial_freqs = c(A_M = 0.10, A_U = 0.40, a_M = 0.10, a_U = 0.40))
    accum_pred     = Theory_accumulative_migration(h_epi = h, m = m, start_E_A = start_E_A, env_B = env_B, env_A = env_A, env_sd = env_sd, n_pop = n_per_pop, t_max = n_gens, drift_C = 10)
    
    # Plot setup
    plot(sim_df$Generation, sim_df$epi_A.mean, type = "n", ylim = c(0, 1),xlab = "Generation", ylab = "Epigenetic mean (pop A)", main = bquote(h[epi] == .(h) ~ ",  m" == .(m)))
    
    # Simulation
    polygon(c(sim_df$Generation, rev(sim_df$Generation)),c(sim_df$epi_A.upr, rev(sim_df$epi_A.lwr)), col = adjustcolor("steelblue", 0.3), border = NA)
    lines(sim_df$Generation, sim_df$epi_A.mean, col = "steelblue", lwd = 2)
    
    # Theory curves
    lines(1:n_gens, simple_pred,       col = "firebrick2",   lty = 1, lwd = 1.5)  # simple
    lines(1:n_gens, migration_pred,    col = "darkorange",   lty = 1, lwd = 1.5)  # migration
    lines(1:n_gens, full_pred$Mean_E,  col = "forestgreen",  lty = 1, lwd = 1.5)  # drift
    lines(1:n_gens, WP_pred$mean_epi,  col = "navy", lty = 1, lwd = 1.5)          # WP model
    lines(1:n_gens, accum_pred$Mean_E, col = "black", lty = 1, lwd = 1.5)         # accumulation
    
    # Legend
    #legend("bottomright", legend = c("Simulation mean", "95% CI", "Theory: simple", "Theory: migration", "Theory: full", "Theory: WP", "Theory: accumulative"), col = c("steelblue", adjustcolor("steelblue", 0.3), "firebrick2", "darkorange", "forestgreen", "navy", "black"), lty = c(1, 1, 1, 1, 1, 1), lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1), bty = "n", cex = 0.75)
  }
}
par(mfrow = c(1, 1))  # Reset



