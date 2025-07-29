data {
  int<lower=0> N_A;               // Number of clusters for Arm A
  int<lower=0> N_B;               // Number of clusters for Arm B
  array[N_A] int<lower=2> cluster_size_A; // Cluster sizes for Arm A (varying)
  array[N_B] int<lower=2> cluster_size_B; // Cluster sizes for Arm B (varying)
  int<lower=2> cluster_size_max;   // Maximum cluster size across all clusters
  array[N_A] vector[cluster_size_max] Y_A; // Outcomes for group A (padded)
  array[N_B] vector[cluster_size_max] Y_B; // Outcomes for group B (padded)
  matrix[cluster_size_max, cluster_size_max] Sigma; // Covariance matrix
}

parameters {
  real mu_A;                    // Mean for group A
  real mu_B;                    // Mean for group B
}

model {
  // Priors
  mu_A ~ normal(0, 1);
  mu_B ~ normal(0, 1);

  // Likelihood (Multivariate Normal with varying sizes)
  for (i in 1:N_A) {
    Y_A[i][1:cluster_size_A[i]] ~ multi_normal(rep_vector(mu_A, cluster_size_A[i]), Sigma[1:cluster_size_A[i], 1:cluster_size_A[i]]);
  }

  for (i in 1:N_B) {
    Y_B[i][1:cluster_size_B[i]] ~ multi_normal(rep_vector(mu_B, cluster_size_B[i]), Sigma[1:cluster_size_B[i], 1:cluster_size_B[i]]);
  }
}
