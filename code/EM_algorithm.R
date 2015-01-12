library(MCMCpack)

initMu <- function(num_features, k){
  Vrnorm <- Vectorize(rnorm)
  
  mu <- Vrnorm(rep(num_features, k))
  mu
}

initSigma <- function(num_features, k){
  Vrinvgamma <- Vectorize(rinvgamma)
  
  sigma <- Vrinvgamma(rep(num_features, k), shape = 1)
  sigma
}

initLambda <- function(num_samples, k){
  #This should be a Dirichlet at some point... but it's late...
  Vrunif <- Vectorize(runif)
  
  lambda <- matrix(nrow=k, ncol=num_samples)
  lambda[1,] <- Vrunif(num_samples)
  
  for(i in 2:(k-1)){
    remaining <- 1 - colSums(lambda, na.rm=T)
    lambda[i,] <- runif(num_samples) * remaining
  }
  
  lambda[k,] <- 1 - colSums(lambda, na.rm=T)
  lambda
}

dataLikelihood <- function(d, mu, sigma, lambda){
  
  sample_mu <- mu %*% lambda
  sample_sigma <- sqrt((sigma ^ 2) %*% lambda)
  
  log_likelihood <- dnorm(d, mean=sample_mu, sd=sample_sigma, log=T)
  log_likelihood <- sum(log_likelihood)
  
  log_likelihood
}

estimateLambda <- function(d, mu, lambda){
  
}

estimateMu <- function(d, mu, lambda, k){
  
  #Init
  s <- dim(mu)
  new_mu <- matrix(nrow=s[1], ncol=s[2])
  
  for(i in 1:k){
    
    temp_mu <- mu[,-i]
    temp_lambda <- lambda[-i,]
    curr_lambda <- lambda[i,]
    
    temp_sample_mu <- temp_mu %*% temp_lambda
    data_discrepancy <- (d - temp_sample_mu) 
    new_mu[,i] <- (data_discrepancy %*% curr_lambda) / (t(curr_lambda) %*% curr_lambda)[1]
  }
  new_mu
}

estimateSigma <- function(d, mu, lambda, sigma){
  
  
}

EMtrain <- function(d, k, num_iterations){
  
  num_features <- nrow(d)
  num_samples <- ncol(d)
  
  mu <- initMu(num_features, k)
  sigma <- initSigma(num_features, k)
  lambda <- initLambda(num_samples, k)
  
  likelihood <- numeric(num_iterations)
  
  for(i in 1:num_iterations){
    #Data Likelihood
    likelihood(i) <- dataLikelihood(d, mu, sigma, lambda)
    
    #E Step
    lambda <- estimateLambda(d, mu, lambda)
    
    #M Step
    mu <- estimateMu(d, mu, lambda, k)
    sigma <- estimateSigma(d, mu, lambda, k, sigma)
  }
  
}