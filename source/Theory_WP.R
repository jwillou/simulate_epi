# Webster & Phillips (2024) - Four epiallele model
# Tracks A_M, A_U, a_M, a_U over generations with mutation and epimutation
Theory_WP <- function(
    N = 100,
    generations = 50,
    mu = 1e-5,
    rho = 1e-3,
    initial_freqs = c(A_M = 0.10, A_U = 0.40, a_M = 0.10, a_U = 0.40)
) {
  var_epi <- numeric(generations)  # to store per-generation individual variance
  # Internal mutation helper
  apply_mutation <- function(allele, mu, rho) {
    parts <- strsplit(allele, "_")[[1]]
    base <- parts[1]
    epi <- parts[2]
    
    if (runif(1) < mu) {
      base <- ifelse(base == "A", "a", "A")
    }
    if (runif(1) < rho) {
      epi <- ifelse(epi == "M", "U", "M")
    }
    paste(base, epi, sep = "_")
  }
  
  `%||%` <- function(x, y) if (!is.null(x)) x else y  # default operator
  
  epialleles <- c("A_M", "A_U", "a_M", "a_U")
  
  get_fitness <- function(geno1, geno2) {
    geno_pair <- paste(sort(c(geno1, geno2)), collapse = "/")
    fitness_map <- list(
      "A_M/A_M" = 1.00,
      "A_M/A_U" = 0.975,
      "A_U/A_U" = 0.95,
      "A_M/a_M" = 0.9,
      "A_U/a_M" = 0.875,
      "A_U/a_U" = 0.85,
      "a_M/a_M" = 0.8,
      "a_M/a_U" = 0.775,
      "a_U/a_U" = 0.70
    )
    return(fitness_map[[geno_pair]] %||% 0.75)
  }
  
  init_pool <- sample(names(initial_freqs), 2 * N, replace = TRUE, prob = initial_freqs)
  population <- matrix(init_pool, nrow = N, ncol = 2, byrow = TRUE)
  
  mean_epi <- numeric(generations)
  freqs <- matrix(0, nrow = generations, ncol = 4)
  colnames(freqs) <- epialleles
  
  for (gen in 1:generations) {
    methylated <- rowSums(matrix(population %in% c("A_M", "a_M"), nrow = N))
    methylation_percent <- methylated / 2
    
    mean_epi[gen] <- mean(methylation_percent)
    var_epi[gen]  <- var(methylation_percent)
    
    all_alleles <- c(population[,1], population[,2])
    freqs[gen, ] <- table(factor(all_alleles, levels = epialleles)) / (2 * N)
    
    fitness_vals <- apply(population, 1, function(ind) get_fitness(ind[1], ind[2]))
    fitness_probs <- fitness_vals / sum(fitness_vals)
    
    parent_indices <- sample(1:N, size = N, replace = TRUE, prob = fitness_probs)
    offspring <- matrix(NA, nrow = N, ncol = 2)
    
    for (i in 1:N) {
      p1 <- population[parent_indices[i], ]
      p2 <- population[sample(1:N, 1), ]
      g1 <- apply_mutation(p1[sample(1:2, 1)], mu, rho)
      g2 <- apply_mutation(p2[sample(1:2, 1)], mu, rho)
      offspring[i, ] <- c(g1, g2)
    }
    
    population <- offspring
  }
  
  result <- as.data.frame(freqs)
  result$Generation <- 1:generations
  result$mean_epi <- mean_epi
  result$var_epi <- var_epi
  return(result)
  
}


# Plot mean epigenetic value
#plot(wp_exact$Generation, wp_exact$mean_epi, type = "l", col = "black", lwd = 2,
#     ylab = "Mean Epigenetic Value", xlab = "Generation", ylim = c(0, 1),
#     main = "W&P Full Simulation (Exact Diploid)")


