library(MASS)
library(lmerTest)

seed_repeat <- 10

avg_power <- rep(NA, seed_repeat)
avg_bias <- rep(NA, seed_repeat)
avg_mse <- rep(NA, seed_repeat)

for (i in 1:seed_repeat) {
  
  seed <- sample(1:1e6, 1)
  set.seed(seed)
  
  n_iter <- 1000
  
  cluster_size <- 4
  icc <- 0.05
  treatment_effect <- 0.1
  
  alpha <- 0.05
  n_clusters_per_arm <- 50
  sigma_w_sq <- 1
  
  results <- rep(NA, n_iter)
  bias <- rep(NA, n_iter)
  mse <- rep(NA, n_iter)
  
  n_clusters <- 2 * n_clusters_per_arm
  total_n <- n_clusters * cluster_size
  cluster_ids <- rep(1:n_clusters, each = cluster_size)
  treatment <- rep(c(0, 1), each = n_clusters_per_arm * cluster_size)
  
  # Between-cluster and within-cluster variance
  sigma_b_sq <- icc / (1 - icc)
  
  # Construct full covariance matrix (block diagonal)
  Sigma <- matrix(0, nrow = total_n, ncol = total_n)
  for (j in 1:n_clusters) {
    idx <- ((j - 1) * cluster_size + 1):(j * cluster_size)
    Sigma[idx, idx] <- sigma_b_sq + diag(sigma_w_sq, cluster_size)
  }
  
  for (m in 1:n_iter){
    mu <- treatment * treatment_effect
    # Simulate one realization
    y <- as.numeric(mvrnorm(n = 1, mu = mu, Sigma = Sigma))
    
    data <- data.frame(
      y = y,
      treatment = factor(treatment),
      cluster = factor(cluster_ids)
    )
    
    # Fit mixed model
    
    model <- lmer(y ~ treatment + (1 | cluster), data = data)
    pval <- summary(model)$coefficients["treatment1", "Pr(>|t|)"]
    est_trt <- summary(model)$coefficients["treatment1", "Estimate"]
    
    bias[m] <- est_trt - treatment_effect
    mse[m] <- (est_trt - treatment_effect)^2
    results[m] <- ifelse(pval < alpha, 1, 0)
    if (m %% 100 == 0) {
      print(paste("Iteration:", m))
    }
  }
  print(paste("Iteration:", i))
  avg_power[i] <- mean(results)
  avg_bias[i] <- mean(bias)
  avg_mse[i] <-mean(mse)
}

mean(avg_power)
mean(avg_bias)
mean(avg_mse)

print(paste0("test_result_", treatment_effect,"_",icc,"_",cluster_size,"_",seed))
write.csv(results, file = paste0("test_result_", treatment_effect,"_",icc,"_",cluster_size,"_",seed))
write.csv(bias, file = paste0("bias_", treatment_effect,"_",icc,"_",cluster_size,"_",seed))
write.csv(mse, file = paste0("mse_", treatment_effect,"_",icc,"_",cluster_size,"_",seed))