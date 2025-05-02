setwd("/Users/jannawilloughby/Desktop/")

# Load the saved RDS
sim_data <- readRDS("output/simulation_output_20250428_100849.rds")

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

#Plot mean simulated epigenetic values with 95% CI shading
h_vals_comp = c(0.01, 0.25, 0.75)
m_vals_comp = m_vals #need to update every time
par(mfrow = c(length(h_vals_comp), length(m_vals_comp)), mar = c(4, 4, 2, 1))
for (h in h_vals_comp) {
  for (m_val in m_vals_comp) {
    df <- subset(agg_df, h_epi == h & m == m_val)
    df <- df[order(df$Generation), ]
    df <- na.omit(df)
    
    plot(df$Generation, df$epi_A.mean, type = "n", ylim = c(0, 1),
         ylab = "Trait Value", xlab = "Generation", main = paste("h_epi =", h, ", m =", m_val))
    
    if (nrow(df) > 1) {
      polygon(c(df$Generation, rev(df$Generation)),
              c(df$epi_A.upr, rev(df$epi_A.lwr)),
              col = rgb(0, 0, 1, 0.2), border = NA)
      polygon(c(df$Generation, rev(df$Generation)),
              c(df$epi_B.upr, rev(df$epi_B.lwr)),
              col = rgb(1, 0, 0, 0.2), border = NA)
      
      lines(df$Generation, df$epi_A.mean, col = "blue", lwd = 2)
      lines(df$Generation, df$epi_B.mean, col = "red", lwd = 2)
    }
    
    legend("topright", legend = c("Epi A", "Epi B"), col = c("blue", "red"), lwd = 2, cex = 0.8)
  }
}

# FST plots with 95% CI shading
par(mfrow = c(length(h_vals_comp), length(m_vals_comp)), mar = c(4, 4, 2, 1))
for (h in h_vals_comp) {
  for (m_val in m_vals_comp) {
    df <- subset(agg_df, h_epi == h & m == m_val)
    df <- df[order(df$Generation), ]
    df <- na.omit(df)
    
    plot(df$Generation, df$epi_FST.mean, type = "n", ylim = c(0, 1),
         ylab = "FST", xlab = "Generation", main = paste("h_epi =", h, ", m =", m_val))
    
    if (nrow(df) > 1) {
      polygon(c(df$Generation, rev(df$Generation)),
              c(df$epi_FST.upr, rev(df$epi_FST.lwr)),
              col = rgb(0, 0.5, 0, 0.2), border = NA)
      polygon(c(df$Generation, rev(df$Generation)),
              c(df$gen_FST.upr, rev(df$gen_FST.lwr)),
              col = rgb(0.5, 0, 0.5, 0.2), border = NA)
      
      lines(df$Generation, df$epi_FST.mean, col = "darkgreen", lwd = 2)
      lines(df$Generation, df$gen_FST.mean, col = "purple", lwd = 2, lty = 2)
    }
    
    legend("topright", legend = c("epi-FST", "gen-FST"), col = c("darkgreen", "purple"),
           lty = c(1, 2), lwd = 2, cex = 0.8)
  }
}

# Reset plotting layout
par(mfrow = c(1, 1))

#Final Generation Data
final_gen <- max(all_df$Generation)

heatmap_df <- all_df %>%
  filter(Generation == final_gen) %>%
  group_by(h_epi, m, rep) %>%
  summarise(
    epi_FST = mean(epi_FST, na.rm = TRUE),
    gen_FST = mean(gen_FST, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(h_epi, m) %>%
  summarise(
    epi_FST_mean = mean(epi_FST, na.rm = TRUE),
    gen_FST_mean = mean(gen_FST, na.rm = TRUE),
    .groups = "drop"
  )

# Reshape matrices
epi_matrix <- with(heatmap_df, tapply(epi_FST_mean, list(h_epi, m), mean))
gen_matrix <- with(heatmap_df, tapply(gen_FST_mean, list(h_epi, m), mean))

# Epigenetic vs Genetic FST Difference
heatmap_df$diff_FST <- heatmap_df$epi_FST_mean - heatmap_df$gen_FST_mean
z_matrix <- with(heatmap_df, tapply(diff_FST, list(h_epi, m), mean))

filled.contour(x = as.numeric(colnames(z_matrix)),
               y = as.numeric(rownames(z_matrix)),
               z = t(z_matrix),  # << again transpose needed
               color.palette = colorRampPalette(c("dodgerblue2", "goldenrod1", "firebrick2")),
               xlab = "Migration rate (m)",
               ylab = "Epigenetic heritability (h_epi)",
               main = "Epi-Gen FST Difference (Epigenetic - Genetic)")


