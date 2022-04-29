data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> Q;
  int<lower=0> K;
  matrix[N, P] y_1;
  matrix[N, Q] y_2;
}
parameters {
  vector[K] z;
  matrix[P, K] W1;
  matrix[Q, K] W2;
  matrix[P, P] invPsi1;
  matrix[Q, Q] invPsi2;
}
model {
  z ~ normal(0, 1);
  W1 ~ normal(0, 1);
  W2 ~ normal(0, 1);
  y1 ~ normal(W1 * z, Psi1)
  y2 ~ normal(W2 * z, Psi2)
}
