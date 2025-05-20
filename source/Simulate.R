Simulate <- function(h, mig_AB, mig_BA, selection, env_A, env_B, env_sd, noise_sd, start_E_A, start_E_B) {
  # carry over object name
  h_epi = h
  
  # Initialize population with SNP genotypes
  geno_A = matrix(rep(x=1, times=(n_per_pop * n_loci * 2)), ncol = 2 * n_loci) #rbinom(n_per_pop * n_loci * 2, 1, env_A)
  geno_B = matrix(rep(x=0, times=(n_per_pop * n_loci * 2)), ncol = 2 * n_loci) #rbinom(n_per_pop * n_loci * 2, 1, env_B)
  
  pop = data.frame(ID = 1:(2 * n_per_pop), pop = rep(c("A", "B"), each = n_per_pop), epi = c(rnorm(n_per_pop, start_E_A, noise_sd), rnorm(n_per_pop, start_E_B, noise_sd)))
  geno = rbind(geno_A, geno_B)
  mean_epi_A = mean_epi_B = numeric(n_gens)
  epi_FST    = gen_FST    = numeric(n_gens)
  
  for (gen in 1:n_gens) {
    # Epigenetic FST
    mean_epi_A[gen] = mean(pop$epi[pop$pop == "A"])
    mean_epi_B[gen] = mean(pop$epi[pop$pop == "B"])
    epi_FST[gen]    = CalcFstEpi(pop$epi, pop$pop)
    
    # Genetic FST
    geno_summed  = sapply(seq(1, ncol(geno), by = 2), function(i) {geno[, i] + geno[, i + 1]})
    geno_summed  = as.data.frame(geno_summed)
    wc_input     = data.frame(pop = ifelse(pop$pop == "A", 1, 2), geno_summed)
    gen_FST[gen] = CalcFstGen(geno_summed, wc_input$pop)
    
    # Reproduction with selection or random mating
    offspring <- list()
    for (pop_group in c("A", "B")) {
      group_idx <- which(pop$pop == pop_group)
      group_pop <- pop[group_idx, ]
      group_geno <- geno[group_idx, ]
      
      # Environmental target for fitness
      env_target <- ifelse(pop_group == "A", env_A, env_B)
      
      if (selection != 0) {
        mismatch = 1 - abs(group_pop$epi - env_target)
        parent_weights = exp(mismatch * selection)
        p1_indices = sample(1:nrow(group_pop), n_per_pop, replace = TRUE, prob = parent_weights)
        p2_indices = sample(1:nrow(group_pop), n_per_pop, replace = TRUE, prob = parent_weights)
      } else {
        p1_indices = sample(1:nrow(group_pop), n_per_pop, replace = TRUE)
        p2_indices = sample(1:nrow(group_pop), n_per_pop, replace = TRUE)
      }
      
      while (any(p1_indices == p2_indices)) {
        same = which(p1_indices == p2_indices)
        p2_indices[same] = sample(1:nrow(group_pop), length(same), replace = TRUE)
      }
      
      # Epigenetic inheritance
      epi_p1 = group_pop$epi[p1_indices]
      epi_p2 = group_pop$epi[p2_indices]
      env_vals = rnorm(n_per_pop, env_target, env_sd)
      epi_new = h_epi * ((epi_p1 + epi_p2) / 2) + (1 - h_epi) * env_vals + rnorm(n_per_pop, 0, noise_sd)
      epi_new = pmin(pmax(epi_new, 0), 1)
      
      # Genotype inheritance
      geno_offspring = matrix(0, nrow = n_per_pop, ncol = ncol(geno))
      for (l in seq(1, ncol(geno), by = 2)) {
        alleles_p1 = ifelse(runif(n_per_pop) < 0.5, group_geno[p1_indices, l], group_geno[p1_indices, l+1])
        alleles_p2 = ifelse(runif(n_per_pop) < 0.5, group_geno[p2_indices, l], group_geno[p2_indices, l+1])
        geno_offspring[, l] = alleles_p1
        geno_offspring[, l+1] = alleles_p2
      }
      
      # Store offspring for this group
      offspring[[pop_group]] <- list(
        epi = epi_new,
        geno = geno_offspring,
        pop = rep(pop_group, n_per_pop)
      )
    }
    
    # Merge A and B offspring
    pop <- data.frame(
      ID = 1:(2 * n_per_pop),
      pop = c(offspring$A$pop, offspring$B$pop),
      epi = c(offspring$A$epi, offspring$B$epi)
    )
    geno <- rbind(offspring$A$geno, offspring$B$geno)
    
    # Migration
    previous_pop <- pop$pop
    is_A = previous_pop == "A"
    is_B = previous_pop == "B"
    new_pop = previous_pop
    migrate_A = runif(sum(is_A)) < mig_AB
    new_pop[which(is_A)[migrate_A]] = "B"
    migrate_B = runif(sum(is_B)) < mig_BA
    new_pop[which(is_B)[migrate_B]] = "A"
    
    # Identify individuals who moved from A into B, "reset" such that B epi is always constant
    t = length(pop$epi[pop$pop=="A" & pop$pop != new_pop])
    if(t>0){
      pop$epi[pop$pop=="A" & pop$pop != new_pop] = rnorm(t, start_E_B, noise_sd)
    }
    
    # Update population identities
    pop$pop = new_pop
    
  }
  data.frame(Generation = 1:n_gens, epi_mean_A = mean_epi_A, epi_mean_B = mean_epi_B, epi_FST = epi_FST, gen_FST = gen_FST)
}
