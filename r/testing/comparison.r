library(gglasso)		# install.package("gglasso")
library(sparseGAM)		# install.packages("sparseGAM")

source("./funcs.r")


# ----------------------------------------
# Test data 
# ----------------------------------------
n <- 200
p <- 1000
gsize <- 10
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(-8, gsize), rep(0, p - gsize))
y <- X %*% b + rnorm(n, 5, .001)



# ----------------------------------------
# Fit methods
# ----------------------------------------
glfit <- gglasso::gglasso(X, y, groups, intercept=FALSE)
ssgl <- SSGL::SSGL(y, X, .1, 300, groups)
ssgl <- sparseGAM::SSGL(y, X, X, groups, lambda0=100)
gsvb <- gsvb::gsvb.fit(y, X, groups, intercept=T, sigma=1e-1)

matplot(t(glfit$beta), type="l")
plot(ssgl$beta)
plot(gsvb$beta)


# ----------------------------------------
# Coefficient est.
# ----------------------------------------
l2(c(5, b) - gsvb$beta_hat)
l2(c(5, b) - c(ssgl$beta0, ssgl$beta))
head(ssgl$beta)
head(f$m * f$g)
head(b)


# ----------------------------------------
# Timing
# ----------------------------------------
microbenchmark::microbenchmark(
    fit(y, X, groups, mu, rep(0.5, p), g, 1, 1, 1.5, 100, 1e-4, T), 
    times=10
)

