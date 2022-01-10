# Mutual Information CCA
library(MASS)

source("./mi_cca_funcs.r")


# ----------------------------------------
# Example 1
# ----------------------------------------
p <- 2; q <- 2; n <- 200
x <- matrix(nrow=n, ncol=p)
y <- matrix(nrow=n, ncol=q)

x[, 1] <- rnorm(n)
x[, 2] <- rnorm(n)

y[, 1] <- x[ , 1]^2 + x[ ,2]^2
y[, 2] <- rnorm(n)

cc <- cca(x, y, cca.objective)

for (i in 1:10) {
    # mi <- cca(x, y, micca.objective_sparse)
    mi <- cca(x, y, micca.objective)
    if (i == 1) maxmi <- mi
    if (mi$val >= maxmi$val) {
	maxmi <- mi
    }
}

plot(maxmi$a %*% t(x), maxmi$b %*% t(y))
plot(cc$a %*% t(x), cc$b %*% t(y))


# ----------------------------------------
# Example 2
# ----------------------------------------
p <- 3; q <- 2; n <- 60
x <- matrix(nrow=n, ncol=p)
y <- matrix(nrow=n, ncol=q)

x[, 1] <- rnorm(n)
x[, 2] <- rchisq(n, 1)
x[, 3] <- rt(n, 5)
y[, 1] <- 2 * x[ , 2] + 3 * x[ , 3]
y[, 2] <- rchisq(n, 2)

cc <- cca(x, y, cca.objective)

for (i in 1:20) {
    mi <- cca(x, y, micca.objective)
    if (i == 1) maxmi <- mi
    if (mi$val >= maxmi$val) {
	maxmi <- mi
    }
}

cc
maxmi


# ----------------------------------------
# Example 3
# ----------------------------------------
p <- 50; q <- 20; n <- 1000
p <- 10; q <- 10; n <- 1000
x <- matrix(nrow=n, ncol=p)
y <- matrix(nrow=n, ncol=q)

x[, 1] <- rnorm(n)
x[, 2] <- rnorm(n)
for (i in 3:p)
    x[ ,i] <- rnorm(n)

y[, 1] <- x[ , 1]^2
y[, 2] <- rnorm(n)
for (i in 3:q)
    y[ ,i] <- rt(n, 4)

cc <- cca(x, y, cca.objective, iter=1e4)

for (i in 1:20) {
    mi <- cca(x, y, micca.objective_sparse, iter=4e2)
    # mi <- cca(x, y, micca.objective)
    if (i == 1) maxmi <- mi
    if (mi$val >= maxmi$val) {
	maxmi <- mi
    }
}

