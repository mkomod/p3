# --- Testing the update eqs. for tau:a and tau:b
library(Rcpp)
Rcpp::sourceCpp("../src/fit.cpp")


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

yty <- sum(y * y)
yx <- t(y) %*% X
xtx <- t(X) %*% X

mu <- rnorm(p)
s <- rgamma(p, 1, 1)
g <- rep(runif(p/5), each=5)


# ----------------------------------------
# Update tau a
# ----------------------------------------
S <- compute_S(yty, yx, xtx, groups, b, rep(0.15, p), !!b, p)

fn <- function(theta, ta0, tb0, S, n) 
{
    ta <- theta[1]
    tb <- theta[2]
    res <- (0.5 * S + tb0 - tb) * (ta / tb) +
	ta * log(tb) - lgamma(ta) +
	(0.5 * n + ta0 - ta) * (log(tb) - digamma(ta))
    return(res)
}

x <- seq(1, 100, 1e-1)
x <- seq(1, 100, 0.5)
y <- x
z <- outer(x, y, Vectorize(function(x, y) update_a_b_obj(x, y, 1, 1, S, n)))
persp(x, y, z, theta=90)
arrayInd(which.min(z), dim(z))
min(z)
update_a_b_obj(43, 43.1, 1, 1, 1e2, n)

optim(c(40, 40), fn=function(x) update_a_b_obj(x[1], x[2], 1e-3, 1e-3, S, n), method="BFGS")
update_a_b(4, 1, 1e-3, 1e-3, S, n)


