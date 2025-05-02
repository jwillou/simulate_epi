Theory_full <- function(h_epi, m, start_E_local, start_E_migrant, env_local, env_migrant, env_sd, n_pop, t_max, drift_C = 1) {
  
  mean_local   = numeric(t_max)
  mean_migrant = numeric(t_max)
  var_local    = numeric(t_max)
  
  mean_local[1]   = start_E_local
  mean_migrant[1] = start_E_migrant
  var_local[1]    = 0.02 
  
  for (t in 2:t_max) {
    # Migration step
    E_after_migration_local   = (1 - m) * mean_local[t-1] + m * mean_migrant[t-1]
    E_after_migration_migrant = (1 - m) * mean_migrant[t-1] + m * mean_local[t-1]
    
    # Deterministic update
    deterministic_local   = h_epi * E_after_migration_local   + (1 - h_epi) * env_local
    deterministic_migrant = h_epi * E_after_migration_migrant + (1 - h_epi) * env_migrant
    
    # Drift based on trait variance (QTL-style sampling error)
    drift_sd_local   = sqrt(drift_C * var_local[t-1] / n_pop)
    drift_sd_migrant = sqrt(drift_C * var_local[t-1] / n_pop)
    
    mean_local[t]   = deterministic_local   + rnorm(1, 0, drift_sd_local)
    mean_migrant[t] = deterministic_migrant + rnorm(1, 0, drift_sd_migrant)
    
    # Update variance (from Lynch & Walsh)
    var_local[t] = (h_epi^2) * var_local[t-1] + (1 - h_epi)^2 * env_sd^2
  }
  
  return(data.frame(Generation = 1:t_max, Mean_E = mean_local, Var_E = var_local))
}

