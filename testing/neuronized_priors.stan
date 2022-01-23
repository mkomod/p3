data {
  int<lower=0> n;
  int<lower=0> p;
  matrix[n, p] X;
  vector[n] y;

  // prior parameters
  real<lower=0> a0;
  real<lower=0> b0;
  real<lower=0> sigma_w;
}
parameters {
  vector<lower=0, upper=1>[p] eta;
  vector[p] w;
  vector[p] a;
  real<lower=0> sigma;
}
transformed parameters {
  vector[p] alpha0;
  vector[p] theta;

  alpha0 = - inv_Phi(eta);
  for (i in 1:p) theta[i] = fmax(0, a[i] - alpha0[i]) * w[i];
}
model {
  eta ~ beta(a0, b0);
  a ~ normal(0, 1);
  w ~ normal(0, sigma_w^2);
  y ~ normal(X * theta, sigma^2);
}
