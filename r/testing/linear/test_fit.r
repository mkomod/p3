library(Rcpp)
sourceCpp("../src/fit.cpp")


# ----------------------------------------
# Test data
# ----------------------------------------
set.seed(1)
n <- 100
p <- 1000
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-4, gsize), rep(8, gsize), rep(0, p - 3 * gsize))
y <- X %*% b + rnorm(n, 0, 1)

# constrained S
f1 <- fit(y, X, groups-1, 1, 1, 200, 1e-3, 1e-3, rnorm(p, 0, 0.2), rep(0.2, p), rep(0.5, p), T,
    T, 1, 500, 150, 1e-3, T)

# unconstrained S
f2 <- fit(y, X, groups-1, 1, 1, 200, 1e-3, 1e-3, rnorm(p, 0, 0.2), rep(0.2, p), rep(0.5, p), F,
    T, 1, 500, 150, 1e-3, T)

f01 <- gsvb.fit(y, X, groups, indp_covaraince=TRUE, niter=15)
f02 <- gsvb.fit(y, X, groups, indp_covaraince=FALSE, niter=15)

# test that w = f2$s is equal to the df/dw
ss <- sapply(f2$S, function(s) sum(diag(matrix(s, nrow=5)))) + 
    sapply(1:200, function(i) sum((f2$m * f2$m)[5*(i-1) + 1:5]))
ww <- as.numeric(f2$s)
mean(((ss)^(-0.5))[groups] - ww)
plot(((ss)^(-0.5))[groups],  ww, xlab="W tilde", ylab="w opt")
lines(c(-5, 10), c(-5, 10))
grid()

plot(f01$beta)
plot(f02$beta)

plot(f01$beta - f02$beta)
plot(b - f02$beta[-1])
sum((b - f01$beta[-1])^2)
sum((b - f02$beta[-1])^2)

plot(f2$m * f2$g) 

sum( (b - (f1$m - f1$g))^2)
sum( (b - (f2$m - f2$g))^2)
plot(f2$g)
plot(f1$g)

plot(f1$s)
sqrt(unlist(lapply(f2$S, function(s) diag(matrix(s, nrow=5)))))
plot(sqrt(unlist(lapply(f2$S, function(s) diag(matrix(s, nrow=5))))))

plot(f1$elbo, type="b")
points(f2$elbo, type="b")
