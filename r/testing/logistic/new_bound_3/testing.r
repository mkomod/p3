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
#
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
g <- f$g[!duplicated(groups)]
plot(sort(g))

gs <- sort(g)
ogs <- order(g)

while(any(cumprod(gs) <= 2.2e-10)) {
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

SS <- as.matrix(S[fs >= 1e-4, ])
rownames(SS) <- NULL
nrow(SS)

m <- f$m
s <- f$s

xm.G <- sapply(unique(groups), function(group) {
    G <- which(groups == group)
    X[ , G] %*% m[G]
})


xs.G <- sapply(unique(groups), function(group) {
    G <- which(groups == group)
    sqrt(X[ , G]^2 %*% s[G]^2)
})


rr <- apply(SS, 1, function(s) {
    j <- ogs[!!s]
    res <- sapply(1:n, function(i) {
	mu <- sum(xm.G[i, j])
	sig <- sqrt(sum(xs.G[i, j]^2))
	new_bound(mu, sig, 101)
    })
    prod(gs^s * (1-gs)^(1-s)) *  sum(res)
})


mid <- which(g <= 0.98 & g >= 0.020)
big <- which(g > 0.98)
SSS <- expand.grid(c(0,1), c(0,1), c(0,1))
gss <- g[mid]

rrr <- apply(SSS, 1, function(s) {
    j <- mid[!!s]
    j <- c(j, big)
    print(s)
    res <- sapply(1:n, function(i) {
	mu <- sum(xm.G[i, j])
	sig <- sqrt(sum(xs.G[i, j]^2))
	new_bound(mu, sig, 101)
    })
    prod(gss^s * (1-gss)^(1-s)) *  sum(res)
})

sum(rrr)



cbind(SS, rr)

sum(rr)

plot(cumsum(sort(rr)))
cbind(S[fs >= 1e-3, ], fs[fs >= 1e-3])

rep(c(0, 1), 14)

