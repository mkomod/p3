# ----------------------------------------
# Test opt S
# ----------------------------------------
set.seed(1)
n <- 10
p <- 5

n <- 10
p <- 20

X <- matrix(rnorm(n * p), nrow=n)
mu <- rnorm(p)

f <- function(S, xtx, mu, X) 
{
    P * sum(exp(X %*% mu + 0.5 * 
    - 0.5 * log(det(S)) + 
    (sum(diag(S) + mu^2))^0.5
}

xtx <- t(X) %*% X

ll <- seq(0.01, 2, by=0.01)
ll <- seq(0.01, 2, by=0.01)

res <- sapply(ll, function(l) {
    S <- solve(xtx + l * diag(rep(1, p)))
    f(S, xtx, mu)
})

plot(ll, res, type="l", lwd=4)


# ----------------------------------------
# x' S x == X'X S [1, 1]?
# ----------------------------------------
n <- 100
p <- 10

X <- matrix(rnorm(n * p), nrow=n)
A <- matrix(rnorm(n * p), nrow=n)
S <- t(A) %*% A

microbenchmark::microbenchmark(
    apply(X, 1, function(x) {
	t(x) %*% S %*% x
    }),
    diag(X %*% S %*% t(X))
)


sum(apply(X, 1, function(x) { t(x) %*% S %*% x }))
diag(X %*% S %*% t(X))
