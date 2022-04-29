# Group sparse Gaussian
n <- 800000
mx1 <- 1;  sx1 <- 4
mx2 <- -3; sx2 <- 0.5
x1 <- rnorm(n, mx1, sx1)
x2 <- rnorm(n, mx2, sx2)

my1 <- 4;  sy1 <- 2
my2 <- 8; sy2 <- 1
gyp <- 0.25
gy <- sample(c(T, F), n, replace=T, prob=c(gyp, 1-gyp))
y1 <- ifelse(gy, rnorm(n, my1, sy1), 0)
y2 <- ifelse(gy, rnorm(n, my2, sy2), 0)

mz1 <- 2;  sz1 <- 1
mz2 <- -2; sz2 <- 0.2
gyz <- 0.8
gz <- sample(c(T, F), n, replace=T, prob=c(gyz, 1-gyz))
z1 <- ifelse(gz, rnorm(n, mz1, sz1), 0)
z2 <- ifelse(gz, rnorm(n, mz2, sz2), 0)

X <- matrix(c(x1 ,x2, y1, y2, z1, z2), ncol=6)

# coviariance matrix
round(cov(X), 2)

# E[b b']
mX <- apply(X, 2, mean)
round(cov(X)+ mX %*% t(mX), 2)


# l2 norm of guassian random vectors
n <- 100000
s <- 2
xs <- matrix(c(rnorm(n, 0, s), rnorm(n, 0, s), rnorm(n, 0, s)), ncol=3)
hist(apply(xs, 1, function(x) sqrt(sum(x^2))))
mean(apply(xs, 1, function(x) sqrt(sum(x^2))))
hist(apply(xs, 1, function(x) sum(x^2)))
p <- 3
2^((1-2)/2) * s * p * gamma((p + 1)/2) / gamma((p + 2)/2)


# quadratic form
s <- c(1, 0.5, 1, 2, 1)^2
mu <- c(0, 4, 1, 5, 1)
b <- diag(1/sqrt(s)) %*% mu
b <- as.numeric(b)
lambda <- s
ssd <- sqrt(s)

xs <- matrix(c(rnorm(n, mu[1], ssd[1]), 
	       rnorm(n, mu[2], ssd[2]),
	       rnorm(n, mu[3], ssd[3]),
	       rnorm(n, mu[4], ssd[4]),
	       rnorm(n, mu[5], ssd[5])), ncol=5)
mean(apply(xs, 1, function(x) sqrt(sum(x^2))))
sqrt(sum(s + mu^2))

U <- matrix(rnorm(n * 3), ncol=3)
Q <- apply(U, 1, function(u) sum(lambda * (u + b)^2))
mean(Q)
mean(sqrt(Q))


# v bound
rx <- matrix(nrow=10, ncol=100)
ax <- matrix(nrow=10, ncol=100)
n <- 50000

for (i in 1:10) {
    for (j in 1:100) {
	mu <- sample(-20:20 * 1/2, i, replace=T) 
	s <- sample(1:20 * 1/5, i, replace=T)^2
	xs <- matrix(rnorm(n * i, mu, sqrt(s)), ncol=i) 
	rx[i, j] <- mean(apply(xs, 1, function(x) sqrt(sum(x^2))))
	ax[i, j] <- sqrt(sum(s + mu^2))
    }
}

plot(apply(rx, 1, mean), ylim=c(0, 25))
points(apply(ax, 1, mean))
plot(rx)
points(ax)



rnorm(50, mu, sqrt(s))

# testing linearlity of inner prod

y <- 1:10
X <- matrix(rnorm(10*5), nrow=10)
b <- rnorm(5)

t(y) %*% X %*% b

sum(sapply(1:5, function(i) t(y) %*% X[ , i] %*% b[i]))

b %*% t(b)
# S <- b %*% t(b)
B <- t(X) %*% X
sum(diag(t(X[, 3:5]) %*% X[, 3:5] %*% S[3:5, 3:5])) + sum(diag(B)[1:2] * diag(S)[1:2])
S <- b %*% t(b)
sum(diag(B %*% S))
t(b) %*% t(X) %*% X %*% b


# sqrt(x)
curve(sqrt, to=1.5, ylim=c(0, 1))
curve(1/4 + x, add=T)
