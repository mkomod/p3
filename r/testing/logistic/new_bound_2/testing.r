# ----------------------------------------
# Multiplication in R
# ----------------------------------------
a <- runif(10)
B <- matrix(runif(100), nrow=10, ncol=10)

a * B == apply(B, 1, function(b) b * a)	 # rows
a * B == apply(B, 2, function(b) b * a)  # cols

# vector is multiplied by column if the same length

a <- runif(10)
B <- matrix(runif(120), nrow=12, ncol=10)

a * B == apply(B, 1, function(b) b * a)	 # rows
a * B == apply(B, 2, function(b) b * a)  # cols

# vector same lenght as row will not raise error if multiplied



# ----------------------------------------
# Testing compute_S
# ----------------------------------------
S <- c()
for (i in 1:n) {
    temp <- 1
    for (group in unique(groups)) {
	G <- which(groups == group)
	temp <- temp * (
	    (1 - g[G][1]) + 
	    g[G][1] * exp(sum(l[i] * X[i, G] * m[G] + 0.5 * 
			      l[i]^2 * X[i, G]^2 * s[G]^2))
	)
    }
    S <- c(S, temp)
}

# ----------------------------------------
# Testing the bound
# ----------------------------------------
# - E log s(-x) <= 
#	E lx + log( E e^(-lx) + E e^(x(1-l)) )
# s(x) = (1 + e^(-x))^-1

s <- function(x) 1/(1 +  exp(-x))

nmgf <- function(m, s, l) exp(l * m + 0.5 * l^2 * s^2)

expec <- function(x, mu, sig) -log(s(-x)) * dnorm(x, mean=mu, sd=sig)

jen <- function(mu, sig) log(1 + nmgf(mu, sig, 1))

jaak <- function(l, mu, sig) {
    -log(s(l)) + (mu + l)/2 + (s(l) - 0.5)/(2*l) * ((sig^2 + mu^2) - l^2)
}

nb0 <- function(mu, sig) 0.5 * (mu + sqrt(2 + (sig^2 + mu^2)))

nb1 <- function(l, mu, sig) {
    l * mu + log(nmgf(mu, sig, -l) + nmgf(mu, sig, 1-l))
}


nb2 <- function(tt, mu, sig) {
    a <- 1 - tt
    b <- 1 + tt

    tt * (sig*sqrt(2/pi)*exp(-mu^2/(2*sig^2)) + mu*(1-2*pnorm(-mu/sig))) +
    log(
    exp(pnorm(mu/sig + sig*-tt, log.p=T)  + -tt*mu + 0.5*tt^2*sig^2) +
    exp(pnorm(-mu/sig + sig*-tt, log.p=T) - -tt*mu + 0.5*tt^2*sig^2) +
    exp(pnorm(mu/sig + sig*a, log.p=T)  + a*mu + 0.5*a^2*sig^2) +
    exp(pnorm(-mu/sig - sig*b, log.p=T) + b*mu + 0.5*b^2*sig^2)
    )
}


nb3 <- function(mu, sig) {
    2 * exp(-mu + sig^2/2) * pnorm(mu/sig - sig) +
    0.5 * sig * sqrt(2/pi) * exp(-mu^2/(2 * sig^2)) +
    mu * pnorm(mu / sig)
}


mu <- 0
sig <- 5

integrate(expec, -50, 50, subdivisions=2e4, mu=mu, sig=sig)
jen(mu, sig)
optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val
nb0(mu, sig)
optim(1, nb1, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val
optim(1, nb2, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val
nb3(mu, sig)

# ----------------------------------------
# Comparison over sig
# ----------------------------------------
mu <- 0
sigs <- seq(1, 10, length.out=40)

res <- sapply(sigs, function(sig) 
{
    c(
      integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value,
      # jen(mu, sig),
      # optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
      # nb0(mu, sig),
      # optim(1, nb1, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
      # optim(1, nb2, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
      nb3(mu, sig)
    )
})


matplot(sigs, t(res), type="l", lwd=3, lty=1)
res

# ----------------------------------------
# Testing NB3
# ----------------------------------------
mu <- 1
sig <- 2

integrate(function(x) (2 * log(1+ exp(-x)) + x) * dnorm(x, mu, sig), 0, Inf)
integrate(function(x) (2 * exp(-x) + x) * dnorm(x, mu, sig), 0, Inf)

integrate(function(x) x * dnorm(x, mean=mu, sd=sig), 0, Inf, subdivisions=1e5)
0.5 * sig * sqrt(2/pi) * exp(-mu^2/(2 * sig^2)) +
mu * pnorm(mu / sig)

integrate(function(x) 2*exp(-x) * dnorm(x, mu, sig), 0, Inf)
2 * exp(-mu + sig^2/2) * pnorm(mu/sig - sig)


# ----------------------------------------
# Simulation
# ----------------------------------------
g <- c(0.5, 0.2, 1, 0.1, 0.9)
m <- -2:2
sig <- rep(0.2, 5)

X <- matrix(rnorm(5), ncol=5)
b <- matrix(rnorm(5000, m, sig) * (runif(5000) <= g), byrow=T, ncol=5)

xb <- X %*% t(b)

xb <- as.numeric(xb)
mean(2 * exp(-xb) * (xb >= 0)) + mean(xb * (xb >= 0))

# ----------------------------------------
# 
# ----------------------------------------
x1 <- rnorm(100000, 5, 1)
x2 <- rnorm(100000, 0, 1)

cov(exp(-x1), x1 >= -x2)

mu <- -11
sig <- 2
integrate(function(x) log(1 + exp(x + mu)) * dnorm(x, 0, sig), -100, 100)
integrate(function(x) log(1+ exp(x)) * dnorm(x, mu, sig), -10, 10)

# this is part (a)
integrate(function(x) (x + mu) * dnorm(x, 0, sig), -mu, Inf)$val
sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * (pnorm(mu/sig))

# (b)
integrate(function(x) exp(mu + x) * dnorm(x, 0, sig), -Inf, -mu)$val
exp(mu + 0.5*sig^2) * pnorm(-mu/sig - sig)

# (c)
integrate(function(x) exp(-(mu + x)) * dnorm(x, 0, sig), -mu, Inf)$val
exp(-mu + 0.5*sig^2) * pnorm(mu/sig - sig)


nb3 <- function(mu, sig) {
    sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * pnorm(mu/sig) +
    exp(mu + 0.5*sig^2) * pnorm(-mu/sig - sig) +
    exp(-mu + 0.5*sig^2) * pnorm(mu/sig - sig)
}


integrate(function(x) log(1 + exp(x)) * dnorm(x, mu, sig), -50, 50)

sigs <- seq(0.5, 10, length.out=40)
pdf(file="~/p1/notes/gsvb/figures/comp_new_bound.pdf", width=9, height=4)
layout(t(1:3))
for (mu in c(-5, 0, 5)) {
    res <- sapply(sigs, function(sig) 
    {
	c(
	  integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value,
	  optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
	  nb3(mu, sig)
	)
    })
    rr <- apply(res[2:3, ], 1,  function(r) log10(r - res[1, ]))
    matplot(sigs, rr, type="l", lwd=3, lty=1, col=c(2,3), xlab=expression(sigma), 
	    ylab="log10 error", main=sprintf("Mean: %d", mu))
    if (mu == -5) {
	legend("topleft", legend=c("Ours", "JJ"), col=c(3,2), lwd=3, lty=c(1))
    }
}
dev.off()

# ----------------------------------------
#
# ----------------------------------------
mus <- seq(-5, 5, length.out=40)
sig <- 0.5
res <- sapply(mus, function(mu) 
{
    c(
      integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value,
      optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
      nb3(mu, sig)
    )
})
layout(1)
rr <- apply(res[2:3, ], 1,  function(r) log10(r - res[1, ]))
matplot(mus, rr, type="l", lwd=3, lty=1, col=c(2,3), xlab=expression(mu), 
	ylab="log10 error", main=sprintf("Sig: %.2f", sig))

G <- which(!!b)
X <- X[ , G]
m <- m[G]
s <- s[G]


jaak <- function(l, mu, sig) {
    -log(sigmoid(l)) + (mu + l)/2 + (sigmoid(l) - 0.5)/(2*l) * ((sig^2 + mu^2) - l^2)
}


tvals <- sapply(1:n, function(i) {
    mu <- sum((X[i, ] + 4) * m)
    sig <- sqrt(sum((X[i, ]^2 + 4) * s^2))
    integrate(function(x) log(1 + exp(x + mu)) * dnorm(x, 0, sig), -100, 100)$val
})


jvals <- sapply(1:n, function(i) {
    mu <- sum((X[i, ] + 4) * m)
    sig <- sqrt(sum((X[i, ] + 4)^2 * s^2))
    optim(1, jaak, mu=mu, sig=sig)$val
})

mus <- (X[, ] + 4) %*% m
nvals <- nb3(mus, sqrt((X[, ] + 4)^2 %*% s^2))

sigmoid <- function(l) 1/(1+exp(-l))


plot(nb, jvals)

sum(nb[nb <= jvals])
sum(nb[nb > jvals])

sum(jvals[nb <= jvals])
sum(jvals[nb > jvals])

plot(nb - tvals, col=1 + (abs(mus) > 1.5))
plot(jvals - tvals, col=1 + (abs(mus) > 1.5))
plot(jvals - nb, col=1 + (abs(mus) > 1.5))

mus <- X[ , Gs] %*% m[Gs]
sigs <- sqrt(X[ , Gs]^2 %*% s[Gs]^2)

jvals <- sapply(1:n, function(i) {
    mu <- mus[i]
    sig <- sigs[i]
    optim(1, jaak, mu=mu, sig=sig)$val
})

nvals <- nb3(mus, sigs)

plot(jvals - nvals)


mu = 0
sig <- 0.7
integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value
optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val
nb3(mu, sig)



