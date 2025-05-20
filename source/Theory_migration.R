Theory_migration <- function(h_epi, m, start_E_local, start_E_migrant, env_local, env_migrant, t_max) {
  E_local      = numeric(t_max)
  E_migrant    = numeric(t_max)
  E_local[1]   = start_E_local
  E_migrant[1] = start_E_migrant
  for (t in 2:t_max) {
    mixed_E_local   = (1 - m) * E_local[t-1] + m * E_migrant[t-1]
    mixed_E_migrant = (1 - m) * E_migrant[t-1] + m * E_local[t-1]
    E_local[t]      = h_epi * mixed_E_local + (1 - h_epi) * env_local
    E_migrant[t]    = h_epi * mixed_E_migrant + (1 - h_epi) * env_migrant
  }
  return(E_local)
}