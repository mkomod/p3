# Mutual Information CCA
library(MASS)


# ----------------------------------------
# Example 1
# ----------------------------------------
p <- 2; q <- 2; n <- 60
x <- matrix(nrow=n, ncol=p)
y <- matrix(nrow=n, ncol=q)

x[, 1] <- rnorm(n)
x[, 2] <- rnorm(n)

y[, 1] <- x[ , 1]^2
y[, 2] <- rnorm(n)

cc <- cca(x, y, cca.objective)

for (i in 1:20) {
    mi <- cca(x, y, micca.objective)
    if (i == 1) maxmi <- mi
    if (mi$val >= maxmi$val) {
	maxmi <- mi
    }
}

maxmi


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

