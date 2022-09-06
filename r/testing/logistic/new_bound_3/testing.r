# ----------------------------------------
# Testing bound #3
# ----------------------------------------

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

jaak <- function(l, mu, sig) {
    -log(sigmoid(l)) + (mu + l)/2 + (sigmoid(l) - 0.5)/(2*l) * ((sig^2 + mu^2) - l^2)
}


tvals <- sapply(1:n, function(i) {
    mu <- sum((X[i, ] + 4) * m)
    sig <- sqrt(sum((X[i, ]^2 + 4) * s^2))
    integrate(function(x) log(1 + exp(x + mu)) * dnorm(x, 0, sig), -100, 100)$val
})

sigmoid <- function(l) 1/(1+exp(-l))

jvals <- sapply(1:n, function(i) 
{
    mu <- sum((X[i, ] + 4) * m)
    sig <- sqrt(sum((X[i, ] + 4)^2 * s^2))
    optim(1, jaak, mu=mu, sig=sig)$val
})



# ----------------------------------------
# bound + variaitonal param
# ----------------------------------------
expec <- function(x, mu, sig) -log(sigmoid(-x)) * dnorm(x, mean=mu, sd=sig)

jaak <- function(l, mu, sig) {
    -log(sigmoid(l)) + (mu + l)/2 + (sigmoid(l) - 0.5)/(2*l) * ((sig^2 + mu^2) - l^2)
}


nb3 <- function(mu, sig) {
    sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * pnorm(mu/sig) +
    exp(mu + 0.5*sig^2) * pnorm(-mu/sig - sig) +
    exp(-mu + 0.5*sig^2) * pnorm(mu/sig - sig)
}


mu <- 0
sig <- 0.5
p <- 1:11
optnb <- function(mu, sig) 
{
    integrate(function(x) (x + mu) * dnorm(x, 0, sig), -mu , Inf)$val +
    integrate(Vectorize(function(x) sum((-1)^(p-1) * exp(p * -(x + mu)) / p) * dnorm(x, 0, sig), "x"), -mu, Inf)$val +
    integrate(Vectorize(function(x) sum((-1)^(p-1) * exp(p * (x + mu)) / p) * dnorm(x, 0, sig), "x"), -Inf, -mu)$val 
}

p <- 1:51
optnb(mu, sig)
nb3(mu, sig)

integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value
optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val
nb3(mu, sig)

xs <- seq(-5, 0, by=.05)
plot(xs, log(1+exp(xs)), type="l")
lines(exp(xs), col=2)

plot(xs, log(1+exp(xs)) - exp(xs))

yx <- exp(xs)
res <- yx - (yx^2)/2 + (yx^3)/3 - yx^4/4 + yx^5/5 - yx^6/6 + yx^7/7 - yx^8/8 + yx^9/9 - yx^10/10 + yx^11/11

p <- 1:101
res <- sapply(xs, function(x) sum((-1)^(p-1) * exp(p * x) / p))


plot(log(1+exp(xs)) - res)


mus <- seq(-5, 5, length.out=40)
sig <- 0.5

res <- sapply(mus, function(mu) 
{
    c(
      integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value,
      optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
      nb3(mu, sig, 9)
    )
})

rr <- apply(res[2:3, ], 1,  function(r) log10(r - res[1, ]))
matplot(mus, rr, type="l", lwd=3, lty=1, col=c(2,3), xlab=expression(mu), 
	ylab="log10 error", main=sprintf("Sig: %.2f", sig))

mus <- seq(-5, 5, length.out=40)
sig <- 5

res <- sapply(mus, function(mu) 
{
    c(
      integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value,
      optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
      nb3(mu, sig, 9)
    )
})

rr <- apply(res[2:3, ], 1,  function(r) log10(r - res[1, ]))
matplot(mus, rr, type="l", lwd=3, lty=1, col=c(2,3), xlab=expression(mu), 
	ylab="log10 error", main=sprintf("Sig: %.2f", sig))


matplot(mus, t(res), type="l", lwd=3, lty=1, col=c(2,3), xlab=expression(mu), 
	ylab="log10 error", main=sprintf("Sig: %.2f", sig))


mu <- 0
sig <- 0.5
nb3 <- function(mu, sig, p=21) {
    sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * pnorm(mu/sig) +
    sum(sapply(1:p, function(l) {
	(-1)^(l-1) / l * exp(mu*l + 0.5*l^2*sig^2 + 
	    pnorm(-mu/sig - l*sig, log.p=T)) +
	(-1)^(l-1) / l * exp(-mu*l + 0.5*l^2*sig^2 + 
	    pnorm(mu/sig - l*sig, log.p=T))
    }))
}

mu <- -5:5
sig <- rep(1, length(mu))

m <- mu[1]
ss <- sig[1]
integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value
optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val
nb3(mu, sig, p=31)

integrate(expec, -50, 50, subdivisions=5e5, mu=m, sig=ss)$value
nb3(m, ss, p=31)

nb3 <- function(mu, sig, p=21) {
    sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * pnorm(mu/sig) +
    sum(sapply(1:p, function(l) {
	(-1)^(l-1) / l * exp(mu*l + 0.5*l^2*sig^2 + 
	    pnorm(-mu/sig - l*sig, log.p=T)) +
	(-1)^(l-1) / l * exp(-mu*l + 0.5*l^2*sig^2 + 
	    pnorm(mu/sig - l*sig, log.p=T))
    }))
}


nb3 <- function(mu, sig, p=21) {
    l <- 1:p
    sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * pnorm(mu/sig) +
    sapply(1:length(mu), function(i) {
	sum(
	(-1)^(l-1) / l * exp(mu[i]*l + 0.5*l^2*sig[i]^2 + 
	    pnorm(-mu[i]/sig[i] - l*sig[i], log.p=T)) +
	(-1)^(l-1) / l * exp(-mu[i]*l + 0.5*l^2*sig[i]^2 + 
	    pnorm(mu[i]/sig[i] - l*sig[i], log.p=T))
	)
    })
}

nb3(mu, sig)
