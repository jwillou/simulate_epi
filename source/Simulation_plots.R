setwd("/Users/jannawilloughby/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/simulate_epi/")
#setwd("/Users/jrw0107/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/simulate_epi/")
directory = getwd()
outdir    = paste(directory,"/output/",sep="")                    #directory to save model output  
source(paste(directory, "/source/FunctionSourcer.R", sep = ''))   #source functions and set source directory

# Load the saved RDS
sim_data = readRDS("../output/simulation_output_20250520_075838.rds")

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

# Pick the combinations to display
h_vals_comp = 0.01#c(0.01, 0.5, 0.9)
m_vals_comp = c(0.01, 0.1, 0.2)

#Trait-mean trajectories with 95 % CI
par(mfrow = c(1,length(h_vals_comp)), mar = c(4, 4, 2, 1))
colors = c("darkolivegreen", "steelblue", "goldenrod") #c("blue", "forestgreen", "orange", "red", "purple")
names(colors) <- as.character(m_vals_comp)

for (h in h_vals_comp) {
  plot(NA, xlim = c(0, max(agg_df$Generation, na.rm = TRUE)), ylim = c(0, 1),  xlab = "Generation", ylab = "Trait value",
       main = sprintf("h = %.2f", h))
  abline(h=0.6, lty=2, col="grey70")
  # Track whether Pop B has been plotted
  pop_B_plotted = FALSE
  for (m_val in m_vals_comp) {
    df <- agg_df %>%
      dplyr::filter(h_epi == h, m_AB == m_val) %>%
      dplyr::arrange(Generation) %>%
      tidyr::drop_na()
    if (nrow(df) > 1) {
      # Only plot Pop B once
      if (!pop_B_plotted) {
        polygon(c(df$Generation, rev(df$Generation)), c(df$epi_B.upr, rev(df$epi_B.lwr)), col = adjustcolor("gray60", alpha.f = 0.3), border = NA)
        lines(df$Generation, df$epi_B.mean, col = "black", lwd = 2, lty = 1)
        pop_B_plotted <- TRUE
      }
      # Plot Pop A polygon (CI) and mean line 
      polygon(c(df$Generation, rev(df$Generation)), c(df$epi_A.upr, rev(df$epi_A.lwr)), col = adjustcolor(colors[as.character(m_val)], alpha.f = 0.2), border = NA)
      lines(df$Generation, df$epi_A.mean, col = colors[as.character(m_val)], lwd = 2)
    }
  }
  #legend("bottomright", legend = c("Pop B (sim)", paste("Pop A, m =", m_vals_comp)),  col = c("black", colors[as.character(m_vals_comp)]), lty = c(2, rep(1, length(m_vals_comp))), lwd = 2, cex = 0.8)
}



# epi-FST vs gen-FST trajectories 

h_vals_comp = 0.01#c(0.01, 0.5, 0.9)
m_vals_comp = c(0.01, 0.1, 0.2)
par(mfrow = c(length(h_vals_comp), length(m_vals_comp)), mar   = c(4, 4, 2, 1))

for (h in h_vals_comp) {
  for (m_val in m_vals_comp) {
    
    df <- agg_df %>%
      dplyr::filter(h_epi == h, m_AB == m_val) %>%
      dplyr::arrange(Generation) %>%
      tidyr::drop_na()
    
    plot(df$Generation, df$epi_FST.mean, type = "n", ylim = c(0, 1), xlab = "Generation", ylab = "FST", main = sprintf("h = %.2f, m = %.3f", h, m_val))
    
    if (nrow(df) > 1) {
      polygon(c(df$Generation, rev(df$Generation)), c(df$epi_FST.upr, rev(df$epi_FST.lwr)), col = rgb(0, 0.5, 0, 0.2), border = NA)
      polygon(c(df$Generation, rev(df$Generation)), c(df$gen_FST.upr, rev(df$gen_FST.lwr)), col = rgb(0.5, 0, 0.5, 0.2), border = NA)
      lines(df$Generation, df$epi_FST.mean, col = "darkgreen", lwd = 2)
      lines(df$Generation, df$gen_FST.mean, col = "purple",    lwd = 2, lty = 2)
    }
    #legend("topright", legend = c("epi-FST", "gen-FST"), col    = c("darkgreen", "purple"), lty    = c(1, 2), lwd = 2, cex = 0.8)
  }
}
par(mfrow = c(1, 1))

# Heat-map of epi-gen FST difference at final generation

final_gen <- max(all_df$Generation)

heatmap_df <- all_df %>%
  dplyr::filter(Generation == final_gen) %>%
  dplyr::group_by(h_epi, m_AB, rep) %>%
  dplyr::summarise(epi_FST = mean(epi_FST, na.rm = TRUE),
                   gen_FST = mean(gen_FST, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::group_by(h_epi, m_AB) %>%
  dplyr::summarise(diff_FST = mean(epi_FST - gen_FST, na.rm = TRUE),
                   .groups = "drop")

z_mat <- with(heatmap_df, tapply(diff_FST, list(h_epi, m_AB), mean))

filled.contour(x = as.numeric(colnames(z_mat)), y = as.numeric(rownames(z_mat)), z = t(z_mat),
               color.palette = colorRampPalette(c("dodgerblue2", "goldenrod1", "firebrick2")),
               xlab = "Migration rate (m_AB)", ylab = "Epigenetic heritability (h_epi)", main = "Mean (epi-FST – gen-FST) at final generation")



