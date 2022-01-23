library(rstan)


# Neuronized priors
wgt <- 1
act <- Vectorize(\(x) max(0, x))

x <- act(rnorm(1e5, 0, 1))
x <- act(rnorm(1e5, 0, 1)) * rnorm(1e5, 0, 1)

hist(x, freq=F, breaks=1000, xlim=c(-3, 3))
hist(x, freq=F, breaks=100)


# Example
n <- 200
p <- 1000
b <- c(1, 1, rep(0, p - 2))
x <- matrix(rnorm(n * p), ncol=p)
y <- x %*% b + rnorm(n, 0, 0.5)

# sdata <- list(n=n, p=p, X=x, y=as.numeric(y), a0=1, b0=1, sigma_w=0.1)
# s <- rstan::stan(file="./neuronized_priors.stan", data=sdata, chains=1)
# plot(s, pars="theta")

sdata <- list(n=n, p=p, X=x, y=as.numeric(y), a0=1, b0=1, sigma_w=0.1)
m <- rstan::stan_model(file="./neuronized_priors.stan")
svb <- rstan::vb(m, data=sdata, iter=2e4, tol_rel_obj=1e-8, algorithm="meanfield")
plot(svb, pars="theta")

