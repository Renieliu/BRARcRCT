library("rstan")
library("blockrand")
library(ggplot2)
library(ggtext)
library(tidyr)
library(rstan)
library(MASS)

# Simulation Parameters
n_cluster <- 100
n_each <- 4
n_repeats <- 500
initial <- 4
iter <- 1000

# True parameters
effectsize <- 0.5
ICC <- 0.15
L <- 0.02
U <- 0.98

true_mu_A <- 1
true_mu_B <- true_mu_A+effectsize
sigma2W <- 1

sigma2B <- ICC / (1 - ICC)
mu_A_vec <- rep(true_mu_A, n_each)
mu_B_vec <- rep(true_mu_B, n_each)

# Storage
allocations2 <- rep(0, n_cluster)
outcomes <- matrix(rep(0, n_cluster*n_each), nrow = n_cluster, ncol = n_each)
y_A <- matrix(nrow = 0, ncol = n_each)
y_B <- matrix(nrow = 0, ncol = n_each)
allocation_plot2 <- matrix(rep(0, n_cluster*n_repeats), nrow = n_cluster, ncol = n_repeats)
test_result <- rep(0, iter)
bias <- matrix(rep(0, iter*n_repeats), nrow = iter, ncol = n_repeats)
mse <- matrix(rep(0, iter*n_repeats), nrow = iter, ncol = n_repeats)

# Covariance matrix for MVN
Sigma <- matrix(sigma2B, nrow = n_each, ncol = n_each)
diag(Sigma) <- sigma2W + sigma2B  # Compound symmetry structure

# Generate block randomized assignment of block size 2 for first 10 participant
randomization <- blockrand(n = initial, block.sizes = 1)
allocations2[1:initial] <- as.character(randomization$treatment)

# Initial allocations (burn-in)
for (i in 1:initial) {
  if (allocations2[i] == "A") {
    ynew <- mvrnorm(n = 1, mu = mu_A_vec, Sigma = Sigma)
    y_A <- rbind(y_A, ynew)
    outcomes[i,] <- ynew
  } else {
    ynew <- mvrnorm(n = 1, mu = mu_B_vec, Sigma = Sigma)
    y_B <- rbind(y_B, ynew)
    outcomes[i,] <- ynew
  }
}

for (m in 1:iter){
  set.seed(m)
  y_A_temp <- y_A 
  y_B_temp <- y_B
  
  for (i in (initial + 1):n_cluster) {
    
    bias_temp <- rep(0, n_repeats)
    mse_temp <- rep(0, n_repeats)
    
    # Stan data format
    stan_data <- list(
      N_A = nrow(y_A_temp),
      N_B = nrow(y_B_temp),
      n_each = n_each,
      Y_A = y_A_temp,
      Y_B = y_B_temp,
      Sigma = Sigma
    )
    
    # Compile and run MCMC
    fit <- stan(file = "BRAR_cRCT.stan", data = stan_data, 
                iter = 2000, chains = 4, refresh = 0)
    # Extract posterior samples
    posterior_samples <- rstan::extract(fit)
    
    for (j in 1:n_repeats) {
      sampled_mu_A <- sample(posterior_samples$mu_A, 1)
      sampled_mu_B <- sample(posterior_samples$mu_B, 1)
      allocation_plot2[i,j] <- ifelse(sampled_mu_A > sampled_mu_B, "A", "B")
      bias_temp[j] <- (sampled_mu_B-sampled_mu_A) - (true_mu_B-true_mu_A)
      mse_temp[j] <- ((sampled_mu_B-sampled_mu_A) - (true_mu_B-true_mu_A))^2
    }
    
    # Final assignment by sampling from the generated allocation distribution
    allocations2[i] <- sample(allocation_plot2[i,], 1)
    
    # Generate new outcome according to allocation and Append to the correct arm
    if (allocations2[i] == "A") {
      ynew <- mvrnorm(n = 1, mu = mu_A_vec, Sigma = Sigma)
      y_A_temp <- rbind(y_A_temp, ynew)
      outcomes[i,] <- ynew
    } else {
      ynew <- mvrnorm(n = 1, mu = mu_B_vec, Sigma = Sigma)
      y_B_temp <- rbind(y_B_temp, ynew)
      outcomes[i,] <- ynew
    }
    
    # print(paste("Iteration:", i))
    prop_B <- sum(allocation_plot2[i,]=="B")/n_repeats
    
    if (prop_B > U) {
      test_result[m] <- 1
      bias[m,] <- bias_temp
      mse[m,] <- mse_temp
      break
    } else if (prop_B < L) {
      test_result[m] <- 0
      bias[m,] <- bias_temp
      mse[m,] <- mse_temp
      break
    } else if (i==n_cluster){
      test_result[m] <- 0
      bias[m,] <- bias_temp
      mse[m,] <- mse_temp
      break
    }
    
  }
  print(paste("Iteration:", m))
}

sum(test_result)/length(test_result)
mean(bias)
mean(mse)
print(paste0("test_result_", effectsize,"_",ICC,"_","(",L,"-",U,")_",n_each))
write.csv(test_result, file = paste0("test_result_", effectsize,"_",ICC,"_","(",L,"-",U,")_",n_each,"_cluster",n_cluster))
write.csv(bias, file = paste0("bias_", effectsize,"_",ICC,"_","(",L,"-",U,")_",n_each,"_cluster",n_cluster))
write.csv(mse, file = paste0("mse_", effectsize,"_",ICC,"_","(",L,"-",U,")_",n_each,"_cluster",n_cluster))