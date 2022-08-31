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
plot(lls, sapply(lls, jaak), type="l")
lines(lls, sapply(lls, nb), type="l")
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

pnorm(0, 3, 2, lower.tail=F)
pnorm((0 - 3)/2, lower.tail=F)

sig <- 2
mu <- 4
tt <- 1
pnorm(0, 2*sig^2*tt - mu, sig, lower.tail=T)
pnorm(mu/sig - 2*sig*tt, 0, 1, lower.tail=T)

xs <- rnorm(100000, 0.2, 4)
mean()

mu <- -4
sig <- 2
tt <- 2
pnorm(mu/sig - 2*sig*tt, 0, 1, lower.tail=T) * exp(2 * tt * (mu + tt * sig^2)) +
pnorm(mu/sig, lower.tail=F)

integrate(function(x) exp((x - abs(x)) * tt) * dnorm(x, mu, sig), -Inf, Inf, 
	  subdivisions=1e4)

