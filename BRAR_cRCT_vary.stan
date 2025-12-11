data {
  int<lower=0> N_A;               // Number of clusters for Arm A
  int<lower=0> N_B;               // Number of clusters for Arm B
  array[N_A] int<lower=2> cluster_size_A; // Cluster sizes for Arm A (varying)
  array[N_B] int<lower=2> cluster_size_B; // Cluster sizes for Arm B (varying)
  int<lower=2> cluster_size_max;   // Maximum cluster size across all clusters
  array[N_A] vector[cluster_size_max] Y_A; // Outcomes for group A (padded)
  array[N_B] vector[cluster_size_max] Y_B; // Outcomes for group B (padded)
  
  // Hyperparameters for variance priors
  real<lower=0> A_sigma_W;        // scale for Half-Cauchy prior on sigma_W
  // For ICC prior: Beta(rho_alpha, rho_beta); set =1,1 for Uniform(0,1)
  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;
}

parameters {
  real mu_A;                    // Mean for group A
  real mu_B;                    // Mean for group B
  real<lower=0> sigma_W;              // Within-cluster SD
  real<lower=0, upper=1> rho;         // ICC
}

transformed parameters {
  real<lower=0> sigma2_W = square(sigma_W);
  real<lower=0> sigma2_B = rho / (1 - rho) * sigma2_W;

  // Common covariance matrix up to cluster_size_max
  matrix[cluster_size_max, cluster_size_max] Sigma;

  // Base: all off-diagonals are sigma2_B
  for (i in 1:cluster_size_max) {
    for (j in 1:cluster_size_max) {
      Sigma[i, j] = sigma2_B;
    }
  }
  // Add within-cluster variance on diagonal
  for (i in 1:cluster_size_max) {
    Sigma[i, i] = Sigma[i, i] + sigma2_W;
  }
}

model {
  // Priors
  mu_A ~ normal(0, 1);
  mu_B ~ normal(0, 1);

  // Half-Cauchy(0, A_sigma_W) on sigma_W (truncated at 0 via <lower=0>)
  sigma_W ~ cauchy(0, A_sigma_W);

  // ICC prior: Beta(rho_alpha, rho_beta)
  rho ~ beta(rho_alpha, rho_beta); 
  // If you want strict Uniform(0,1), just set rho_alpha = rho_beta = 1 from R.

  // Likelihood (Multivariate Normal with varying sizes)
  for (i in 1:N_A) {
    Y_A[i][1:cluster_size_A[i]] ~ multi_normal(
      rep_vector(mu_A, cluster_size_A[i]),
      Sigma[1:cluster_size_A[i], 1:cluster_size_A[i]]
    );
  }

  for (i in 1:N_B) {
    Y_B[i][1:cluster_size_B[i]] ~ multi_normal(
      rep_vector(mu_B, cluster_size_B[i]),
      Sigma[1:cluster_size_B[i], 1:cluster_size_B[i]]
    );
  }
}
