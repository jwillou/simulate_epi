#setwd("/Users/jannawilloughby/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/")
setwd("/Users/jrw0107/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/")
directory = getwd()
outdir    = paste(directory,"/output/",sep="")                    #directory to save model output  
source(paste(directory, "/source/FunctionSourcer.R", sep = ''))   #source functions and set source directory

#Parameters
set.seed(2112)
n_per_pop = 100
n_gens    = 50
n_reps    = 100
h_vals    = 0.5 #c(0.01, 0.25, 0.5, 0.75, 1)
m_AB_vals = 0.5 #c(0.01, 0.2, 0.5)
m_BA_val  = 0   #does not iterate over this one, should be only 1 value
env_A     = 0.2
env_B     = 0.8
env_sd    = 0
noise_sd  = 0
n_loci    = 50
selection = 0 #0 turns off selection

# Create progress log file
logfile <- "simulation_progress.txt"
write("Starting simulations...\n", file = logfile)

#Run replicates
results_list <- list()
for (h in h_vals) {
  for (m in m_AB_vals) {
    reps = replicate(n_reps, Simulate(h, m, m_BA_val, selection, env_A, env_B, env_sd), simplify = FALSE)
    reps_df = do.call(rbind, lapply(seq_along(reps), function(i) {
      cbind(reps[[i]], rep = i)
    }))
    reps_df$h_epi = h
    reps_df$m_AB = m
    reps_df$m_BA = m_BA_val
    results_list[[paste(h, m, m_BA_val, sep = "_")]] = reps_df
    
    # Log progress
    cat(sprintf("Finished h = %.2f, m_AB = %.2f, m_BA = %.2f at %s\n", h, m, m_BA_val, format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = logfile, append = TRUE)
  }
}

# Combine and aggregate
all_df <- do.call(rbind, results_list)
agg_df <- all_df %>%
  group_by(Generation, h_epi, m) %>%
  summarise(
    epi_A.mean = mean(epi_mean_A, na.rm = TRUE),
    epi_A.lwr  = quantile(epi_mean_A, 0.025, na.rm = TRUE),
    epi_A.upr  = quantile(epi_mean_A, 0.975, na.rm = TRUE),
    
    epi_B.mean = mean(epi_mean_B, na.rm = TRUE),
    epi_B.lwr  = quantile(epi_mean_B, 0.025, na.rm = TRUE),
    epi_B.upr  = quantile(epi_mean_B, 0.975, na.rm = TRUE),
    
    epi_FST.mean = mean(epi_FST, na.rm = TRUE),
    epi_FST.lwr  = quantile(epi_FST, 0.025, na.rm = TRUE),
    epi_FST.upr  = quantile(epi_FST, 0.975, na.rm = TRUE),
    
    gen_FST.mean = mean(gen_FST, na.rm = TRUE),
    gen_FST.lwr  = quantile(gen_FST, 0.025, na.rm = TRUE),
    gen_FST.upr  = quantile(gen_FST, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Save simulated outputs, parameters, log
timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")
saveRDS(list(all_df = all_df, agg_df = agg_df, n_per_pop = n_per_pop, n_gens = n_gens, n_reps = n_reps, h_vals = h_vals, m_vals = m_vals, env_A = env_A, env_B = env_B, env_sd = env_sd, noise_sd = noise_sd, n_loci = n_loci), file = paste0("../output/simulation_output_", timestamp, ".rds"))
file.copy(logfile, paste0("../output/simulation_progress_", timestamp, ".txt"))

