Theory_accumulative_migration <- function(h_epi, m, start_E_A, env_B, env_A, env_sd, n_pop, t_max, drift_C) {
  mu_A  = numeric(t_max)  # mean trait in population A
  q_B   = numeric(t_max)  # proportion of ancestry from B
  var_A = numeric(t_max)  # variance in trait values for A
  
  # Initialize
  mu_A[1]  = start_E_A
  q_B[1]   = 0
  var_A[1] = env_sd^2  # assume equilibrium env-induced variance at t=1
  
  for (t in 2:t_max) {
    # Update B ancestry proportion with heritability retention
    q_B[t] = (1 - m) * (h_epi * q_B[t - 1]) + m
    
    # Inherited component based on migrant ancestry
    inherited_component = (1 - q_B[t]) * mu_A[t - 1] + q_B[t] * env_B
    
    # Drift variance term (Lynch & Walsh-style)
    drift_sd = sqrt(drift_C * var_A[t - 1] / n_pop)
    
    # Update trait mean with stochastic sampling
    mu_A[t] = h_epi * inherited_component + (1 - h_epi) * env_A + rnorm(1, mean = 0, sd = drift_sd)
    
    # Update variance recursively
    var_A[t] = (h_epi^2) * var_A[t - 1] + (1 - h_epi)^2 * env_sd^2
  }
  
  return(data.frame(Generation = 1:t_max, Mean_E = mu_A, Var_E = var_A, qB = q_B))
}
