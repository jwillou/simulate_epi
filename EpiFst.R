setwd("/Users/jannawilloughby/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/simulate_epi/")
#setwd("/Users/jrw0107/Google Drive/My Drive/Willoughby lab/projects - active/epi migration/simulate_epi/")
directory = getwd()
outdir    = paste(directory,"/output/",sep="")                    #directory to save model output  
source(paste(directory, "/source/FunctionSourcer.R", sep = ''))   #source functions and set source directory

#Parameters
set.seed(2112)
n_per_pop = 100                          # Match theory n_pop=100
n_gens    = 50                           # Match theory t_max=250
n_reps    = 100
h_vals    = c(0.01, 0.5, 0.75, 0.9, 1)   # Match theory compare_h_vals=c(0.5, 0.75, 0.9) 0.01, 0.5, 0.75, 0.9, 1
m_AB_vals = c(0.01, 0.05, 0.1, 0.2)      # Match theory compare_m_vals=c(0.01, 0.1)  0.01, 0.1, 0.2
m_BA_val  = 0                            # Allow option for migration to be asymmetric, currently set symmetric below
env_A     = 0.5                          # Match theory 0.5
env_B     = 1                            # Match theory 1
env_sd    = 0.1                          # Match theory environmental noise 0.02
start_E_A = 0.2                          # Match theory 0.2
start_E_B = 1                            # Match theory 1
noise_sd  = 0.01                         # No epi noise beyond env_sd=0
n_loci    = 250
selection = 0.2                          # No selection=0

# Create progress log file
logfile <- "simulation_progress.txt"
write("Starting simulations...\n", file = logfile)

#Run replicates
results_list <- list()
for (h in h_vals) {
  for (m in m_AB_vals) {
    m_BA_val = m                   # symmetric migration
    reps = replicate(n_reps, Simulate(h, m, m_BA_val, selection, env_A, env_B, env_sd, noise_sd, start_E_A, start_E_B), simplify = FALSE)
    reps_df = do.call(rbind, lapply(seq_along(reps), \(i) cbind(reps[[i]], rep = i)))
    reps_df$h_epi = h
    reps_df$m_AB  = m
    reps_df$m_BA  = m_BA_val
    results_list[[paste(h, m, m_BA_val, sep = "_")]] = reps_df
    
    cat(sprintf("Finished h = %.2f, m_AB = %.2f, m_BA = %.2f at %s\n",  h, m, m_BA_val, format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = logfile, append = TRUE)
  }
}

# aggregate output into manageble dataframe
all_df <- do.call(rbind, results_list)
agg_df <- all_df %>%
  group_by(Generation, h_epi, m_AB) %>%
  summarise(
    epi_A.mean  = mean(epi_mean_A, na.rm = TRUE),
    epi_A.lwr   = quantile(epi_mean_A, .025, na.rm = TRUE),
    epi_A.upr   = quantile(epi_mean_A, .975, na.rm = TRUE),
    
    epi_B.mean  = mean(epi_mean_B, na.rm = TRUE),
    epi_B.lwr   = quantile(epi_mean_B, .025, na.rm = TRUE),
    epi_B.upr   = quantile(epi_mean_B, .975, na.rm = TRUE),
    
    epi_FST.mean = mean(epi_FST, na.rm = TRUE),
    epi_FST.lwr  = quantile(epi_FST, .025, na.rm = TRUE),
    epi_FST.upr  = quantile(epi_FST, .975, na.rm = TRUE),
    
    gen_FST.mean = mean(gen_FST, na.rm = TRUE),
    gen_FST.lwr  = quantile(gen_FST, .025, na.rm = TRUE),
    gen_FST.upr  = quantile(gen_FST, .975, na.rm = TRUE),
    .groups = "drop"
  )

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")   
saveRDS(list(all_df   = all_df,
             agg_df   = agg_df,
             n_per_pop = n_per_pop,
             n_gens    = n_gens,
             n_reps    = n_reps,
             h_vals    = h_vals,
             m_AB_vals = m_AB_vals,
             m_BA_val  = m_BA_val,
             env_A     = env_A,
             env_B     = env_B,
             env_sd    = env_sd,
             noise_sd  = noise_sd,
             n_loci    = n_loci),
        file = file.path(outdir, paste0("simulation_output_", timestamp, ".rds")))
file.copy(logfile, file.path(outdir, paste0("simulation_output_", timestamp, ".txt")))
file.copy(file.path(directory, "EpiFst.R"), file.path(outdir, paste0("simulation_output_", timestamp, "EpiFst.txt")))




