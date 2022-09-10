# ----------------------------------------
# Testing bound #3
# ----------------------------------------
expec <- function(x, mu, sig) -log(sigmoid(-x)) * dnorm(x, mean=mu, sd=sig)

jaak <- function(l, mu, sig) {
    -log(sigmoid(l)) + (mu + l)/2 + (sigmoid(l) - 0.5)/(2*l) * ((sig^2 + mu^2) - l^2)
}

sigmoid <- function(x) 1/(1+exp(-x))

new_bound <- function(mu, sig, p=21) {
    l <- 1:p
    sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * pnorm(mu/sig) +
    sum(
	(-1)^(l-1) / l * exp(mu*l + 0.5*l^2*sig^2 + 
	    pnorm(-mu/sig - l*sig, log.p=T)) +
	(-1)^(l-1) / l * exp(-mu*l + 0.5*l^2*sig^2 + 
	    pnorm(mu/sig - l*sig, log.p=T))
    )
}


# ----------------------------------------
# Testing comupation
# ----------------------------------------
mu <- 1
sig <- 2

# part (a)
integrate(function(x) (x + mu) * dnorm(x, 0, sig), -mu, Inf)$val
sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * (pnorm(mu/sig))

# part (b)
integrate(function(x) exp(mu + x) * dnorm(x, 0, sig), -Inf, -mu)$val
exp(mu + 0.5*sig^2) * pnorm(-mu/sig - sig)

# part (c)
integrate(function(x) exp(-(mu + x)) * dnorm(x, 0, sig), -mu, Inf)$val
exp(-mu + 0.5*sig^2) * pnorm(mu/sig - sig)


# ----------------------------------------
# Test computation when p > 1
# ----------------------------------------
new_bound_int <- function(mu, sig, p) 
{
    p <- 1:p
    integrate(function(x) (x + mu) * dnorm(x, 0, sig), -mu , Inf)$val +
    integrate(Vectorize(function(x) {
	sum((-1)^(p-1) * exp(p * -(x + mu)) / p) * dnorm(x, 0, sig)
    }, "x"), -mu, Inf)$val +
    integrate(Vectorize(function(x) {
	sum((-1)^(p-1) * exp(p * (x + mu)) / p) * dnorm(x, 0, sig)
    }, "x"), -Inf, -mu)$val 
}

new_bound_int(mu, sig, 21)
new_bound(mu, sig, 21)
integrate(expec, -100, 100, mu=mu, sig=sig)


# ----------------------------------------
# Test compuation C++
# ----------------------------------------
Rcpp::sourceCpp("../new_bound_3/bound.cpp")

tll(1:10, 1:10, 30)
sum(sapply(1:10, function(i) new_bound(i, i, 30*2 - 1)))

# about x2 faster for C++
microbenchmark::microbenchmark(
    tll(1:10, 1:10, 30),
    sum(sapply(1:10, function(i) new_bound(i, i, 30*2 - 1)))
)

tll(1:10, 1:10, 30) ==
sum(sapply(1:10, function(i) new_bound(i, i, 30*2 - 1)))


# ----------------------------------------
# Comparison to other bounds
# ----------------------------------------
library(latex2exp)

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
	    ylab=cat(expression(log[10]), " error"), main=sprintf("Mean: %d", mu))
    if (mu == -5) {
	legend("topleft", legend=c("Ours", "Jaakkola"), col=c(3,2), lwd=3, lty=c(1))
    }
}
dev.off()


mus <- seq(-5, 5, length.out=40)
sig <- 1


pdf(file="~/p1/notes/gsvb/figures/comp_new_bound_sig.pdf", width=9, height=4)
    layout(t(1:3))
    for (sig in c(1, 3, 5)) 
    {
	res <- sapply(mus, function(mu) {
	    c(integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value,
	      optim(1, jaak, method="Brent", lower=-4, upper=4, mu=mu, sig=sig)$val,
	      new_bound(mu, sig, 15))
	})
	rr <- apply(res[2:3, ], 1,  function(r) log10(r - res[1, ]))

	matplot(mus, rr, type="l", lwd=3, lty=1, col=c(2,3), xlab=expression(mu),
	    ylim=c(-9, 0),
	    ylab=TeX("$\\log_{10}$ error"), main=sprintf("std. dev: %d", sig))
	if (sig == 1) {
	    legend("topleft", legend=c("Ours", "Jaakkola"), 
		   col=c(3,2), lwd=3, lty=c(1))
	}
    }
dev.off()


# ----------------------------------------
# Testing number of terms needed
# ----------------------------------------
mus <- seq(-5, 5, by=0.5)
sigs <- seq(0.5, 5, by=0.5)
ms <- expand.grid(mus, sigs)
ps <- c()

for (i in 1:nrow(ms)) {
    mu <- ms[i, 1]
    sig <- ms[i, 2]
    for (p in 1:150) {
	if (
	    abs(integrate(expec, -50, 50, subdivisions=5e5, mu=mu, sig=sig)$value - 
		new_bound(mu, sig, 2*p-1)) <= 1e-2
	) {
	    ps <- c(ps, 2*p - 1)
	    break
	}
    }
}

pdf(file="~/p1/notes/gsvb/figures/opt_l.pdf", width=6, height=4)
    par(mar=c(4, 4, 1, 1), family="Times")
    plot(NULL, ylim=c(0.5, 5), xlim=c(-5, 5),
	 ylab=expression(sigma), xlab=expression(mu), axes=F)
    for (i in seq_along(ps)) points(ms[i, 1], ms[i, 2], pch=as.character((ps[i] + 1)/2),
	col=c("green", "blue", "orange", "darkorange", "red")[(ps[i] + 1)/2])
    abline(v=c(-5.5, mus, 5.5)+0.25, lwd=0.2)
    abline(h=c(0.07, sigs, 4.93)+0.25, lwd=.2)
    axis(1, at = -5:5, labels = -5:5, tick=F)
    axis(2, at = sigs, las = 1, labsels=sigs, tick=F)
dev.off()     


# ----------------------------------------
# Testing with gamma
# ----------------------------------------
set.seed(1)
g <- c(runif(90, 0, 0.1), runif(5, 0.95, 1), runif(5, 0.1, 0.9))
plot(sort(g))

g <- f$g[!duplicated(groups)]
gs <- sort(g)
ogs <- order(g)

while(any(cumprod(gs) <= 2.2e-10)) 
{
    ogs <- ogs[cumprod(gs) <= 2.2e-10]
    gs <- gs[cumprod(gs) <= 2.2e-10]
}

sort(g)

S <- expand.grid(
    c(0,1), c(0,1), c(0,1), c(0,1),
    c(0,1), c(0,1), c(0,1), c(0,1),
    c(0,1), c(0,1), c(0,1), c(0,1),
    c(0,1), c(0,1)
)

fs <- apply(S, 1, function(s) {
    prod(gs^s * (1-gs)^(1-s))
})


xm.G <- sapply(unique(groups), function(group) {
    G <- which(groups == group)
    X[ , G] %*% m[G]
})

xs.G <- sapply(unique(groups), function(group) {
    G <- which(groups == group)
    sqrt(X[ , G]^2 %*% s[G]^2)
})

mid <- which(g <= 0.98 & g >= 0.020)
big <- which(g > 0.98)
tot <- 0
gk <- g[mid]

for (i in 0:(2^length(mid)-1)) 
{
    sk <- as.integer(intToBits(i)[1:3])
    j <- c(mid[!!sk], big)

    res <- sapply(1:n, function(i) 
    {
	mu <- sum(xm.G[i, j])
	sig <- sqrt(sum(xs.G[i, j]^2))
	new_bound(mu, sig, 101)
    })
    
    tot <- tot + prod(gk^sk * (1 - gk)^(1 - sk)) *  sum(res)
}


# ----------------------------------------
# Testing computation for SpSL
# ----------------------------------------
e_ll <- function(X.m, X.s, g, tau, groups)
{
    g <- g[!duplicated(groups)]
    mid <- which(g >= tau & g <= (1-tau))
    big <- which(g > (1-tau))
    gk <- g[mid]
    tot <- 0

    for (i in 0:(2^length(mid)-1)) 
    {
	sk <- as.integer(intToBits(i)[1:length(mid)])
	J <- c(mid[!!sk], big)
	mu <- apply(X.m[ , J], 1, sum)
	sig <- sqrt(apply(X.s[ , J], 1, sum))

	res <- tll(mu, sig, 101)
	tot <- tot + prod(gk^sk * (1 - gk)^(1 - sk)) *  sum(res)
    }

    tot
}

Xm <- sapply(unique(groups), function(group) {
    G <- which(groups == group)
    X[ , G] %*% m[G]
})

Xs <- sapply(unique(groups), function(group) {
    G <- which(groups == group)
    X[ , G]^2 %*% s[G]^2
})

e_ll(Xm, Xs, g, 0.025, groups)


# ----------------------------------------
# Testing optimization
# ----------------------------------------
e_ll <- function(X.m, X.s, ug, tau, l=20)
{
    mid <- which(ug >= tau & ug <= (1-tau))
    big <- which(ug > (1-tau))
    gk <- ug[mid]
    tot <- 0

    if (length(mid) >= 10) print("shittles")
    
    if (length(mid) == 0) {
	if (length(big) == 1) {
	    mu <- X.m[ , big]
	    sig <- sqrt(X.s[ , big])
	} else {
	    mu <- apply(X.m[ , big], 1, sum)
	    sig <- sqrt(apply(X.s[ , big], 1, sum))
	}
	return(tll(mu, sig, l))
    }

    for (i in 0:(2^length(mid)-1)) 
    {
	sk <- as.integer(intToBits(i)[1:length(mid)])
	J <- c(mid[!!sk], big)
	
	if (length(J) == 1) {
	    mu <- X.m[ , J]
	    sig <- sqrt(X.s[ , J])
	} else {
	    mu <- apply(X.m[ , J], 1, sum)
	    sig <- sqrt(apply(X.s[ , J], 1, sum))
	}

	res <- tll(mu, sig, l)
	tot <- tot + prod(gk^sk * (1 - gk)^(1 - sk)) *  sum(res)
    }

    tot
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

opt_mu <- function(m_G, y, X, m, s, g, G, lambda, tau, p) 
{
    m[G] <- m_G 
    J <- union(G, which(g >= tau))

    mu <- X[ , J] %*% m[J]
    sig <- sqrt(X[ , J]^2 %*% s[J]^2)

    sum(nb3(mu, sig, p))  -
    sum(y * (X[ , G] %*% m[G])) +
    lambda * sqrt(sum(s[G]^2) + sum(m_G^2))
}

opt_mu_ell <- function(m_G, y, X, m, s, ug, lambda, group, G,
    X.m, X.s, thresh=0.02, l=20) 
{
    m[G] <- m_G
    ug[group] <- 1

    xm <- X[ , G] %*% m[G]
    X.m[ , group] <- xm

    e_ll(X.m, X.s, ug, thresh, l) -
    sum(y * xm) +
    lambda * sqrt(sum(s[G]^2) + sum(m_G^2))
}


m <- f$m
s <- f$s
g <- f$g
g <- +!!b
ug <- g[!duplicated(groups)]

X.m <- sapply(unique(groups), function(group) {
    G <- which(groups == group)
    X[ , G] %*% m[G]
})

X.s <- sapply(unique(groups), function(group) {
    G <- which(groups == group)
    X[ , G]^2 %*% s[G]^2
})

G <- 1:5
group <- 1

opt_mu(1:5, y, X, m, s, g, G, lambda, 0.02, 41)
opt_mu_ell(1:5, y, X, m, s, ug, lambda, group, G, X.m, X.s)

microbenchmark::microbenchmark(
    optim(m[G], 
	fn=function(mG) opt_mu_ell(mG, y, X, m, s, ug, lambda, group,
		G, X.m, X.s),
	control=list(maxit=30),
	method="BFGS"),
    optim(m[G], 
	fn=function(mG) opt_mu(mG, y, X, m, s, g, G, lambda, 0.02, 41),
	control=list(maxit=30),
	method="BFGS")
)


# ----------------------------------------
# Testing opt_g
# ----------------------------------------
opt_g <- function(y, X, m, s, ug, lambda, group, G, X.m, X.s,
	tau=0.02, l=20)
{
    mk <- length(G)
    Ck <- mk * log(2) + 0.5 * (mk -1) * log(pi) + lgamma(0.5 * (mk + 1))

    ug[group] <- 1
    S1 <- e_ll(X.m, X.s, ug, tau, l)

    ug[group] <- 0
    S0 <- e_ll(X.m, X.s, ug, tau, l)

    res <- 
	log(w / (1 - w)) + 
	0.5 * mk -
	Ck +
	mk * log(lambda) +
	0.5 * sum(log(2 * pi * s[G]^2)) -
	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) - 
	S1 + sum(y * (X[ , G] %*% m[G])) + S0

    sigmoid(res)
}


group <- 4
G <- which(groups == group)
opt_g(y, X, m, s, ug, lambda, group, G, X.m, X.s, tau)

group <- 1
ug[group] <- 1
G <- which(groups == group)
opt_mu(m[G], y, X, m, s, ug, lambda, group, G, X.m, X.s, 0.02, 10)

group <- 1
ug[group] <- 1
e_ll(X.m, X.s, ug, tau, l)


# ----------------------------------------
# Testing derivatives
# ----------------------------------------
Rcpp::sourceCpp("../new_bound_3/bound.cpp")

# testing against finite differences - PASS
mu <- X[ , G] %*% m[G]
sig <- sqrt(X[ , G]^2 %*% s[G]^2)
t1 <- tll(mu, sig, 20)

m.h <- m[G]
m.h[1] <- m.h[1] - 1e-10
mu.h <- X[ , G] %*% m.h[G]
t2 <- tll(mu.h, sig, 20)

(t1 - t2) / 1e-10
dt_dm(X, mu, sig, G-1, 20)[1]

# testing for larger group - PASS
mu <- apply(X.m[ , 1:3], 1, sum)
sig <- sqrt(apply(X.s[ , 1:3], 1, sum))
t1 <- tll(mu, sig, 20)

GG <- 1:15
m.h <- m[GG]
m.h[1] <- m.h[1] - 1e-10
mu.h <- apply(X[ , GG] %*% m.h[GG], 1, sum)
t2 <- tll(mu.h, sig, 20)

(t1 - t2) / 1e-10
dt_dm(X, mu, sig, G-1, 20)

# testing sigma
mu <- X[ , G] %*% m[G]
sig <- sqrt(X[ , G]^2 %*% s[G]^2)
t1 <- tll(mu, sig, 20)

s.h <- s[G]
s.h[1] <- s.h[1] - 1e-10
sig.h <- sqrt(X[ , G]^2 %*% s.h[G]^2)
t2 <- tll(mu, sig.h, 20)

(t1 - t2) / 1e-10
dt_ds(X, s, mu, sig, G-1, 20)[1]

