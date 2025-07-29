data {
  int<lower=0> N_A;          // Number of clusters for Arm A
  int<lower=0> N_B;          // Number of clusters for Arm B
  int<lower=0> n_each;
  matrix[N_A, n_each] Y_A;      // Outcomes matrix for group A
  matrix[N_B, n_each] Y_B;      // Outcomes matrix for group B
  matrix[n_each, n_each] Sigma;  // Covariance matrix
}

parameters {
  real mu_A;                    // Mean for group A
  real mu_B;                    // Mean for group B
}

model {
  // Priors
  mu_A ~ normal(0, 1);
  mu_B ~ normal(0, 1);

  // Likelihood (Multivariate Normal)
  for (i in 1:N_A) {
    Y_A[i,] ~ multi_normal(rep_vector(mu_A, n_each), Sigma);
  }

  for (i in 1:N_B) {
    Y_B[i,] ~ multi_normal(rep_vector(mu_B, n_each), Sigma);
  }
}
