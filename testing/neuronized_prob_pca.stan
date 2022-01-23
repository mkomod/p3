data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0> k;
  matrix[n, p] X;

  // prior parameters
  real a0;
  real<lower=0> sigma_w;
}
parameters {
  matrix[n, k] z;
  real<lower=0> sigma_x;

  matrix[p, k] a;
  matrix[p, k] w;
}
transformed parameters {
  matrix[p, k] B;

  for (i in 1:p) {
    for (j in 1:k) {
      B[i, j] = fmax(0, a[i, j] - a0) * w[i, j];
    }
  }
}
model {
  to_vector(a) ~ normal(0, 1);
  to_vector(w) ~ normal(0, sigma_w);

  to_vector(z) ~ normal(0, 1);
  to_vector(X) ~ normal(to_vector(z * B'), sigma_w^2);
}
