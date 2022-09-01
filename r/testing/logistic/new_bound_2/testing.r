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
nb <- function(l, mu, sig) l * mu + log(nmgf(mu, sig, -l) + nmgf(mu, sig, 1-l))
jaak <- function(l, mu, sig) -log(s(l)) + (mu + l)/2 + (s(l) - 0.5)/(2*l) * ((sig^2 + mu^2) - l^2)

nb.4 <- function(l, k, mu, sig)
{
    k * (sig^2 + mu^2) + l * mu + log(nmgf(mu, sig, -l) + nmgf(mu, sig, 1-l))
}

mu <- 2
sig <- 0.25
x <- rnorm(100000, mean=mu, sd=sig)

integrate(expec, -50, 50, subdivisions=2e4, mu=mu, sig=sig)
(mcs <- mean(-log(s(-x))))
(jen <- log(1 + nmgf(mu, sig, 1)))
(nb1 <- 0.5 * (mu + sqrt(2 + (sig^2 + mu^2))))
optim(1, nb, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val
optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val


lls <- seq(-5, 5, length.out=200)
plot(lls, sapply(lls, jaak, mu=mu, sig=sig), type="l")
lines(lls, sapply(lls, nb, mu=mu, sig=sig), type="l")
abline(h=jen)
abline(h=mcs)

# comparison of the bounds with change in mu
mus <- seq(-10, 10, length.out=40)
sigs <- seq(0.2, 6, length.out=40)

sig <- 0.5
# res <- sapply(sigs, function(sig) 
res <- sapply(mus, function(mu) 
{
    x <- rnorm(20000, mean=mu, sd=sig)
    c(
      integrate(expec, -10, 50, subdivisions=5e4, mu=mu, sig=sig)$value,
      mean(-log(s(-x))),
      log(1 + nmgf(mu, sig, 1)),
      optim(1, nb, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
      optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val
    )
})

ldiff <- log(apply(res, 1, function(r) abs(r - res[1, ])))
ldiff[ , 2] <- smooth.spline(sigs, ldiff[ , 2])$y
matplot(sigs, ldiff, type="l", lwd=3, col=c(0, 4:1), lty=1)
legend("bottomright", legend=c("MCI", "Jen", "NB", "JJ"), col=4:1, lwd=3, lty=1)

matplot(mus, ldiff, type="l", lwd=3, col=c(0, 4:1), lty=1)
legend("bottomright", legend=c("MCI", "Jen", "NB", "JJ"), col=4:1, lwd=3, lty=1)

X <- matrix(rnorm(4e6, 0, 1), nrow=400, ncol=1e3)
b <- c(rep(2, 5), rep(-2, 5), rep(0, 990))
Xb <- X %*% b

mean(Xb)
sd(Xb)

nb2 <- function(l, mu, sig) l * (sig^2 + mu^2) + 
    log(nmgf(mu, sig, -l) + nmgf(mu, sig, 1-l))


xs <- seq(-3, 5, length.out=500)
fx <- exp(xs) * dnorm(xs)
plot(xs, fx, type="l")
mean(xs)


mu <- 1
sig <- 2

x <- rnorm(20000, mean=mu, sd=sig)
integrate(expec, -10, 50, subdivisions=5e4, mu=mu, sig=sig)$value
mean(-log(s(-x)))

nb3 <- function(tt, mu, sig) {
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

nb3.int <- function(tt, mu, sig) {
    integrate(function(x) tt * abs(x) * dnorm(x, mu, sig), -1e2, 1e2,
	      subdivisions=1e5)$val + 
    log(
    integrate(function(x) exp((x - tt * abs(x))) * dnorm(x, mu, sig), -1e2, 1e2,
	      subdivisions=1e5)$val +
    integrate(function(x) exp((-tt * abs(x))) * dnorm(x, mu, sig), -1e2, 1e2,
	  subdivisions=1e5)$val
    )
}


mu <- 1
sig <- 0.5

optim(0.5, nb3, method="Brent", lower=-3, upper=3, mu=mu, sig=sig)$val
optim(0.5, nb3.int, method="Brent", lower=-3, upper=3, mu=mu, sig=sig)$val
optim(1, nb, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val
optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val


# mus <- seq(-10, 10, length.out=40)
mu <- 0
sigs <- seq(1, 100, length.out=40)

sig <- 15
mu <- 0.5
res <- sapply(sigs, function(sig) 
{
    x <- rnorm(20000, mean=mu, sd=sig)
    c(
      integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value,
      # mean(-log(s(-x))),
      # log(1 + nmgf(mu, sig, 1)),
      # optim(1, nb, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
      optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
      nb.5(mu, sig)
      # optim(0.5, nb3, method="Brent", lower=-10, upper=10, mu=mu, sig=sig)$val
      # optim(0.5, nb3.int, method="Brent", lower=-3, upper=3, mu=mu, sig=sig)$val
    )
})



res
rr <- apply(res[-1, ], 1, function(r) r - res[1, ])
matplot(sigs, log(rr), type="l")

matplot(mus, log(abs(rr)), type="l")
matplot(mus, t(res), type="l")
matplot(sigs, t(res), type="l")

lls <- seq(-5, 5, length.out=200)
plot(lls, sapply(lls, jaak, mu=mu, sig=sig), type="l")
lines(lls, sapply(lls, nb, mu=mu, sig=sig), type="l")
lines(lls, sapply(lls, nb3, mu=mu, sig=sig), type="l")
sapply(lls, nb, mu=mu, sig=sig)
sapply(lls, nb3, mu=mu, sig=sig)
abline(h=jen)
abline(h=mcs)

nb.2 <- function(x) 0.5 * (x + sqrt(2 + x^2)) + x*tt

nb.4 <- function(theta, mu, sig)
{
    l <- theta[1]
    k <- theta[2]
    k * (sig^2 + mu^2) + l * mu + log(nmgf(mu, sig, -l) + nmgf(mu, sig, 1-l))
}

optim(c(1, 0), nb.4, mu=mu, sig=sig)

nb.5 <- function(mu, sig) {
    2 * exp(-mu + sig^2/2) * pnorm(mu/sig - sig) +
    0.5 * sig * sqrt(2/pi) * exp(-mu^2/(2 * sig^2)) +
    mu * pnorm(mu / sig)
}

integrate(function(x) (2 * log(1+ exp(-x)) + x) * dnorm(x, mu, sig), 0, Inf)
integrate(function(x) (2 * exp(-x) + x) * dnorm(x, mu, sig), 0, Inf)
nb.5(mu, sig)


optim(3, jaak, mu=mu,sig=sig)$val
mu <- 0
sig <- 6
integrate(expec, -50, 50, subdivisions=2e4, mu=mu, sig=sig)
nb.5(mu, sig)

mu <- 1
sig <- 2
integrate(function(x) x * dnorm(x, mean=mu, sd=sig), 0, Inf, subdivisions=1e5)
0.5 * sig * sqrt(2/pi) * exp(-mu^2/(2 * sig^2)) +
mu * pnorm(mu / sig)

integrate(function(x) 2*exp(-x) * dnorm(x, mu, sig), 0, Inf)
2 * exp(-mu + sig^2/2) * pnorm(mu/sig - sig)



