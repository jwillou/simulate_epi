CalcFstGen <- function(geno_summed, pop_ids) {
  fst_per_locus = apply(geno_summed, 2, function(x) {
    group_means = tapply(x, pop_ids, mean)
    overall_mean = mean(x)
    V_B = sum(tapply(x, pop_ids, function(xi) length(xi) * (mean(xi) - overall_mean)^2)) / length(x)
    V_W = mean(tapply(x, pop_ids, var))
    if ((V_B + V_W) == 0) return(NA)
    return(V_B / (V_B + V_W))
  })
  return(mean(fst_per_locus, na.rm = TRUE))
}