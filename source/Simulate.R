Simulate <- function(h, mig_AB, mig_BA, selection, env_A, env_B, env_sd) {
  h_epi = h
  
  # Initialize population with SNP genotypes
  geno_A <- matrix(rbinom(n_per_pop * n_loci * 2, 1, env_A), ncol = 2 * n_loci)
  geno_B <- matrix(rbinom(n_per_pop * n_loci * 2, 1, env_B), ncol = 2 * n_loci)
  
  pop <- data.frame(
    ID = 1:(2 * n_per_pop),
    pop = rep(c("A", "B"), each = n_per_pop),
    epi = c(rnorm(n_per_pop, env_A, 0.01), rnorm(n_per_pop, env_B, 0.01))
  )
  geno <- rbind(geno_A, geno_B)
  mean_epi_A <- mean_epi_B <- numeric(n_gens)
  epi_FST <- gen_FST <- numeric(n_gens)
  
  for (gen in 1:n_gens) {
    # Epigenetic FST
    mean_epi_A[gen] <- mean(pop$epi[pop$pop == "A"])
    mean_epi_B[gen] <- mean(pop$epi[pop$pop == "B"])
    epi_FST[gen]    <- CalcFst(pop$epi, pop$pop)
    
    # Genetic FST
    geno_summed <- sapply(seq(1, ncol(geno), by = 2), function(i) {
      geno[, i] + geno[, i + 1]
    })
    geno_summed <- as.data.frame(geno_summed)
    wc_input <- data.frame(pop = ifelse(pop$pop == "A", 1, 2), geno_summed)
    wc_result <- wc(wc_input)
    gen_FST[gen] <- wc_result$FST
    
    # Reproduction with selection or random mating
    env_target <- ifelse(pop$pop == "A", env_A, env_B)
    if (selection != 0) {
      mismatch <- 1 - abs(pop$epi - env_target)
      parent_weights <- exp(mismatch * selection)
      p1_indices <- sample(1:nrow(pop), nrow(pop), replace = TRUE, prob = parent_weights)
      p2_indices <- sample(1:nrow(pop), nrow(pop), replace = TRUE, prob = parent_weights)
    } else {
      p1_indices <- sample(1:nrow(pop), nrow(pop), replace = TRUE)
      p2_indices <- sample(1:nrow(pop), nrow(pop), replace = TRUE)
    }
    
    # No selfing
    while (any(p1_indices == p2_indices)) {
      same <- which(p1_indices == p2_indices)
      p2_indices[same] <- sample(1:nrow(pop), length(same), replace = TRUE)
    }
    
    # Epigenetic inheritance
    epi_p1 <- pop$epi[p1_indices]
    epi_p2 <- pop$epi[p2_indices]
    current_pop <- pop$pop[p1_indices]
    env_vals <- ifelse(current_pop == "A",
                       rnorm(nrow(pop), env_A, env_sd),
                       rnorm(nrow(pop), env_B, env_sd))
    epi_new <- h_epi * ((epi_p1 + epi_p2) / 2) + (1 - h_epi) * env_vals + rnorm(nrow(pop), 0, noise_sd)
    epi_new <- pmin(pmax(epi_new, 0), 1)
    
    # Genotype inheritance
    geno_offspring <- matrix(0, nrow = nrow(pop), ncol = ncol(geno))
    for (l in seq(1, ncol(geno), by = 2)) {
      alleles_p1 <- ifelse(runif(nrow(pop)) < 0.5, geno[p1_indices, l], geno[p1_indices, l + 1])
      alleles_p2 <- ifelse(runif(nrow(pop)) < 0.5, geno[p2_indices, l], geno[p2_indices, l + 1])
      geno_offspring[, l] <- alleles_p1
      geno_offspring[, l + 1] <- alleles_p2
    }
    
    # Reset all B-born individuals to match env_B
    born_in_B <- current_pop == "B"
    epi_new[born_in_B] <- rnorm(sum(born_in_B), env_B, env_sd)
    
    # Migration AFTER reproduction
    new_pop <- current_pop
    
    # Identify individuals in A and B
    is_A <- current_pop == "A"
    is_B <- current_pop == "B"
    
    # For A → B migration
    migrate_A <- runif(sum(is_A)) < mig_AB
    new_pop[which(is_A)[migrate_A]] <- "B"
    
    # For B → A migration
    migrate_B <- runif(sum(is_B)) < mig_BA
    new_pop[which(is_B)[migrate_B]] <- "A"
    
    # Identify individuals who moved into B and "reset" such that B epi is always constant
    migrated_to_B <- (new_pop == "B") & (current_pop != "B")
    epi_new[migrated_to_B] <- rnorm(sum(migrated_to_B), env_B, env_sd)
    
    # Check overall population size, correct as needed
    offspring_df <- data.frame(
      epi = epi_new,
      pop = new_pop,
      stringsAsFactors = FALSE
    )
    
    # Maintain fixed size: n_per_pop in A and B
    pop_A_idx <- which(offspring_df$pop == "A")
    pop_B_idx <- which(offspring_df$pop == "B")
    
    # Sample to enforce fixed population sizes
    if (length(pop_A_idx) >= n_per_pop) {
      keep_A <- sample(pop_A_idx, n_per_pop, replace = FALSE)
    } else {
      # Not enough A individuals — sample with replacement to reach 100
      keep_A <- sample(pop_A_idx, n_per_pop, replace = TRUE)
    }
    
    if (length(pop_B_idx) >= n_per_pop) {
      keep_B <- sample(pop_B_idx, n_per_pop, replace = FALSE)
    } else {
      # Not enough B individuals — sample with replacement to reach 100
      keep_B <- sample(pop_B_idx, n_per_pop, replace = TRUE)
    }
    
    # Combine final population
    final_indices <- c(keep_A, keep_B)
    pop <- data.frame(
      ID = 1:(2 * n_per_pop),
      pop = offspring_df$pop[final_indices],
      epi = offspring_df$epi[final_indices]
    )
    geno <- geno_offspring[final_indices, ]
    
    # Final population and genotype update
    pop$epi <- epi_new
    pop$pop <- new_pop
    geno <- geno_offspring
  }
  
  data.frame(Generation = 1:n_gens,
             epi_mean_A = mean_epi_A,
             epi_mean_B = mean_epi_B,
             epi_FST = epi_FST,
             gen_FST = gen_FST)
}
