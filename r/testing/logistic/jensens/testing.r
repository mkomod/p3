library(Rcpp)

Rcpp::sourceCpp("./jensens.cpp")

approxEq <- function(a, b, e=1e-3) all(abs(a - b) < e)

# ----------------------------------------
# Test data
# ----------------------------------------
set.seed(1)
n <- 350
p <- 500
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-1, gsize), rep(1, gsize), rep(0, p - 3 * gsize))
xb <- X %*% b
prob <- 1 / (1 + exp(-xb))
y <- rbinom(n, 1, prob)



# ----------------------------------------
# MVN MGF
# ----------------------------------------
n_mgf <- function(X, m, s)
{
    apply(X, 1, function(x) {
	exp(sum(x * m + 0.5 * x^2 * s^2))
    })
}

approxEq(mvnMGF(X, b, rep(.1, p)), n_mgf(X, b, rep(.1, p)))


test_minus(matrix(1:10, ncol=2), 1:5)
