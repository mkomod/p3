data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> K;
  matrix[N, P] x;

  // prior parameters
  real<lower=0> x_a0;
  real<lower=0> x_b0;
  real<lower=0> W_a0;
  real<lower=0> W_b0;
}
parameters {
  matrix[N, K] z;
  matrix[P, K] W;
  real<lower=0> tau_x;
  vector<lower=0>[K] tau_W;
}
transformed parameters {
  real<lower=0> sigma_x;
  vector<lower=0>[K] sigma_W;
  sigma_x = inv(sqrt(tau_x));
  sigma_W = inv(sqrt(tau_W));
}
model {
  tau_x ~ gamma(x_a0, x_b0);
  tau_W ~ gamma(W_a0, W_b0);
  for (k in 1:K) W[ , k] ~ normal(0, sigma_W[k]) ;

  to_vector(z) ~ normal(0, 1);
  to_vector(x) ~ normal(to_vector(z * W'), sigma_x);
}

