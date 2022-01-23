library(rstan)

# Test data
n <- 100
p <- 5
k <- 3
B <- matrix(c(
	0,  2,  0, 
	0,  0,  2, 
	4,  0,  0, 
	0, -1,  0, 
	1,  0,  0
     ), ncol=k, nrow=p)
z <- matrix(rnorm(n*k), ncol=k)
X <- z %*% t(B)

sdata <- list(n=n, p=p, k=k, X=X, a0=0, sigma_w=0.5)
s <- rstan::stan(file="./neuronized_prob_pca.stan", data=sdata, chains=1)

B <- rstan::extract(s, pars="B")[[1]]
plot(s, pars="B")
round(apply(B, c(2, 3), function(b) mean(b[!!b])), 2)

matplot(B[ , , 1], type="l", ylim=c(-5, 5))
matplot(B[ , , 2], type="l", add=T)
matplot(B[ , , 3], type="l", col=1:5, lty=1, add=T)


m <- rstan::stan_model("./neuronized_prob_pca.stan")
sdata <- list(n=n, p=p, k=k, X=X, a0=0, sigma_w=0.2)
svb <- rstan::vb(m, sdata)
B <- rstan::extract(svb, pars="B")[[1]]
plot(svb, pars="B")
round(apply(B, c(2, 3), function(b) mean(b[!!b])), 1)

matplot(B[ , , 1], type="l")
matplot(B[ , , 2], type="l", add=T)
matplot(B[ , , 3], type="l", col=1:5, lty=1, add=T)

