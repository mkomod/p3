# Misc code for debugging
source("./00-functions.R")

d <- dgp_diag(200, 1000, 5, 10, list(corr=0))
m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
f <- m_gsvb(d)
f <- gsvb::gsvb.fit(d$y, d$X, d$groups, intercept=TRUE, a0=1, b0=200, track_elbo=FALSE)

d <- dgp_block(200, 1000, 5, 3, list(corr=0.6, block_size=50), seed=84)
m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
f <- m_gsvb(d, m_par=m_par)

d <- dgp_block(200, 1000, 5, 3, list(corr=0.6, block_size=50), seed=91)
m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
f <- m_gsvb(d, m_par=m_par)

f <- function(x) 
{
    a <- tryCatch({
	ftime <- system.time(x)
	TRUE
    }, error=function(e) {
	return(FALSE)
    })

    if (!a) {
	return(NA)
    }
    return(ftime)
}

f(stop("123"))
f(2)


d <- dgp_wishart(200, 1000, 5, 3, list(dof=3, wieght=0.9))
d <- dgp_diag(200, 1000, 5, 3, pars=list(corr=0))
m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
f <- m_gsvb(d)
f <- gsvb::gsvb.fit(d$y, d$X, d$groups, intercept=TRUE, a0=1, b0=200, track_elbo=FALSE)


# ----------------------------------------
# Posterior predictive
# ----------------------------------------
# and posterior predictive coverage in 95% int?
d$X <- cbind(1, d$X)
d$groups <- c(1, d$groups+1)

y. <- seq(-20, 20, by=0.05)
x. <- d$X[1, ]
res <- replicate(10000, {
    b <- rnorm(length(f$m), f$m, f$s) * (runif(length(f$g)) <= f$g)[d$groups]
    1 / sqrt(2 * pi) * 
	exp(lgamma(f$tau_a + 0.5) - lgamma(f$tau_a)) *
	exp(
	    f$tau_a * log(f$tau_b) +
	    (f$tau_a + 0.5) * (log(2) - log( (y. - sum(x. * b))^2 + 2 * f$tau_b))
	)
}, simplify="matrix")

sig <- sqrt(f$tau_b / f$tau_a)
ys <- replicate(100, {
    # b <- rnorm(length(f$m), f$m, f$s) * (runif(length(f$g)) <= f$g)[d$groups]
    sum(x.* f$beta_hat) + sig * rt(1, 2*f$tau_a)
}, simplify="matrix")

ys <- sum(x.* f$beta_hat) + sig * rt(10000, 2*f$tau_a)

plot(density(ys))
lines(y., probs , type="l", lwd=4)
lines(density(ystrue) , type="l", lwd=4)
abline(v=d$y[1])



probs <- apply(res, 1, mean)
plot(y., probs , type="l", lwd=4)
abline(v=d$y[1])
matplot(y., res, type="l")





cvrg <- apply(res, 1, quantile, probs=c(0.025, 0.975))

quants <- apply(res, 1, quantile, probs=c(0.005, 0.995))
plot(d$y, d$X %*% f$beta_hat)



mean((quants[1, ] <= d$y) & (d$y <= quants[2, ]))
mean((cvrg[1, ] <= c(0, d$b)) & (c(0, d$b) <= cvrg[2, ]))


mf <- spsl::spsl.group_sparse(d$y, d$X, d$groups)
mf$T <- mf$T[, -(1:500)]


ystrue <- sapply(1:4500, function(i) {
    b <- mf$B[ , i] * mf$Z[ , i][d$groups]
    tau <- mf$T[i]
    rnorm(1, sum(x. * b), sqrt(tau))
})

plot(density(ys), xlim=c(-10, 10))
lines(density(ystrue), xlim=c(-10, 10))
plot(density(ystrue))
lines(density(ys))
abline(v=d$y[1])



# ----------------------------------------
# Coverage?
# ----------------------------------------





