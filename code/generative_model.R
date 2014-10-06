library(ggplot2)

#Generative model
universe <- rnorm(5000, 2, 0.75)
qplot(universe, geom='histogram')

dirichlet_sample <- function(n, samples, c){
  
  #Determining beta dist
  num_samples <- length(samples)
  b <- rbeta(num_samples, 1, c)
  #Init
  p <- numeric(num_samples)
  p[1] <- b[1]
  
  #Stick-breaking process
  p[2:n] <- sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
  
  theta <- sample(samples, size=n, prob=p, replace=TRUE)
  theta
}

tau_dist <- dirichlet_sample(3000, universe, 10)
qplot(tau_dist, geom='histogram')

type_dist <- dirichlet_sample(100, tau_dist, 10)
qplot(type_dist, geom='histogram')

num_samples <- 10

determine_samples <- function(distribution, n){
  #Init
  sample_values <- c()
  
  #Figuring out a weighted average of each expression level's contribution
  exp_levels <- sort(unique(distribution))
  prob_weights <- as.numeric(table(distribution))/ length(distribution)
  exp_weights <- exp_levels * prob_weights
  
  for (i in 1:n){
    #Determining how much of each type contributes to the sample
    sample_type_amount <- rnorm(length(exp_levels), 1, 0.5)
    sample_total <- sample_type_amount %*% exp_weights
    
    #Adding it to the list
    sample_values <- c(sample_values, sample_total)
  }
  sample_values
}

single_gene_samples <- determine_samples(type_dist, 10)
qplot(single_gene_samples, geom='density')

simulate_data <- function(num_genes, num_samples){
  universe <- rnorm(5000, 2, 0.75)
  
  samples <- c()
  for(i in 1:num_genes){
    
    tau_dist <- dirichlet_sample(2000, universe, 10)
    type_dist <- dirichlet_sample(500, tau_dist, 10)
    
    gene_samples <- determine_samples(type_dist, num_samples)
    samples <- rbind(samples, gene_samples)
  }
  rownames(samples) <- 1:length(samples[,1])
  samples
}
