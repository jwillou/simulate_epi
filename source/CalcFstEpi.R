#Function to compute FST from trait data
CalcFstEpi <- function(trait, pop_ids) {
  group_means  = tapply(trait, pop_ids, mean)
  overall_mean = mean(trait)
  V_B          = sum(tapply(trait, pop_ids, function(x) length(x) * (mean(x) - overall_mean)^2)) / length(trait)
  V_W          = mean(tapply(trait, pop_ids, var))
  return(V_B / (V_B + V_W))
}