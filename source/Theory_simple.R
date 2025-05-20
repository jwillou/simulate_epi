Theory_simple <- function(h_epi, start_E, env_target, t_max) {
  E    = numeric(t_max)
  E[1] = start_E
  for (t in 2:t_max) {
    E[t] = h_epi * E[t-1] + (1 - h_epi) * env_target
  }
  return(E)
}