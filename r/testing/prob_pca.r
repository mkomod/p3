z <- matrix(rnorm(2 * 200), ncol=2)
W <- matrix(c(2, 2, 0, 1, -2, 4), ncol=2)
x <- t(W %*% t(z))



s0_data <- list(N=nrow(x), D=ncol(x), K=3, X=x)
s0 <- rstan::stan(model_code=pca, data=s0_data,
	    chains=1, iter=4e3)

s1_data <- list(N=nrow(x), P=ncol(x), K=3, x=x, 
	    x_a0=1, x_b0=1, W_a0=1e-3, W_b0=1e-3)
s1 <- rstan::stan(file="./prob_pca.stan", data=s1_data,
	    chains=1, iter=4e3)

plot(s0, pars="W")
plot(s1, pars="W")

s_data <- list(N=nrow(x), D=ncol(x), K=2, X=x)

m <- rstan::stan_model(file="./prob_pca.stan")
s.vb <- rstan::vb(m, data=s_data, algorithm="meanfield")

s.params <- extract(s)
plot(s, pars="W")

SVD <- svd(cov(x))
SVD$u
SVD$u / min(abs(SVD$u))
