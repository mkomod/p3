# Misc code for debugging
source("./00-functions.R")

d <- dgp_diag(250, 1000, 10, 10, list(model="gaussian", corr=0))
d <- dgp_diag(250, 1000, 10, 10, list(model="gaussian", corr=0))
d <- dgp_diag(300, 1000, 5, 3, list(model="binomial", corr=0))
d <- dgp_diag(250, 1000, 5, 3, list(model="poisson", corr=0))

plot(d$y)

m_par <- list(family="gaussian", lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3, diag_covariance=TRUE)
m_par <- list(family="binomial-jensens", lambda=1, a0=1, b0=200, diag_covariance=FALSE)
m_par <- list(family="binomial-jaakkola", lambda=1, a0=1, b0=200, diag_covariance=FALSE)
m_par <- list(family="binomial-jaakkola", lambda=1, a0=1, b0=200, diag_covariance=TRUE)
m_par <- list(family="binomial-refined", lambda=1, a0=1, b0=200, diag_covariance=FALSE)
m_par <- list(family="poisson", lambda=1, a0=1, b0=200, diag_covariance=FALSE)

f <- m_gsvb(d, m_par)
f <- m_spsl(d)
f <- m_ssgl(d)

# ----------------------------------------
# m_run
# ----------------------------------------
setting_parameters <- list(n=250, p=1e3, g=5, s=3, dgp=dgp_wishart,
	pars=list(model="gaussian", dof=3, weight=0.9), runs=10)

m_par <- list(family="gaussian", lambda=1, a0=1, b0=200, 
    a_t=1e-3, b_t=1e-3, diag_covariance=TRUE)
m_run(m_gsvb, m_par, setting_parameters, 2)

m_par <- list(family="gaussian", l0=100, l1=1, a0=1, b0=200)
m_run(m_ssgl, m_par, setting_parameters, 2)




# ----------------------------------------
# tryCatch
# ----------------------------------------
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


# ----------------------------------------
# VB posterior predictive
# ----------------------------------------
d <- dgp_diag(200, 1000, 5, 3, pars=list(corr=0))
f <- gsvb::gsvb.fit(d$y, d$X, d$groups, intercept=TRUE, a0=1, b0=200, track_elbo=TRUE)

gsvb.predict <- function(fit, groups, newdata, samples=1e4, quantiles=c(0.025, 0.975), return_samples=FALSE) 
{
    if (fit$parameters$intercept) {
	newdata <- cbind(1, newdata)
	groups <- c(1, groups + 1)
    }

    M <- length(fit$g)
    n <- nrow(newdata)
    sigma <- sqrt(f$tau_b / f$tau_a)

    y.star <- replicate(samples, 
    {
	grp <- (runif(M) <= fit$g)[groups]
	mu <- rnorm(sum(grp), fit$mu[grp], fit$s[grp])
	newdata[ , grp] %*% mu + sigma * rt(n, 2 * f$tau_a)
    }, simplify="matrix")

    y.star <- matrix(y.star, nrow=n)

    res <- list(
	mean=apply(y.star, 1, mean),
	quantiles=apply(y.star, 1, quantile, probs=quantiles)
    )

    if (return_samples) {
	res <- c(res, list(samples=y.star))
    }

    return(res)
}



pred <- gsvb.predict(f, d$groups, d$X, 1e4, return_samples=F)
# coverage
mean(
    (pred$quantiles[1, ] <= d$y) & (pred$quantiles[2, ] >= d$y)
)

pred <- gsvb.predict(f, d$groups, matrix(d$X[1, ], nrow=1), 1e3, return_samples=T)
plot(density(pred$samples))
abline(v=d$y[1])

# ---------------------------------------- 
# Posterior predictive SpSL
# ---------------------------------------- 
d <- dgp_diag(200, 1000, 5, 3, pars=list(corr=0))
mf <- spsl::spsl.group_sparse(d$y, d$X, d$groups, intercept=TRUE, mcmc_samples=2e4)
mf$T <- mf$T[ , -(1:2000)]

spsl.predict <- function(fit, groups, newdata, quantiles=c(0.025, 0.975), return_samples=FALSE) 
{
    if (fit$settings$intercept) {
	groups <- c(1, groups+1)
	newdata <- cbind(1, newdata)
    }

    n <- nrow(newdata)
    tau <- sqrt(fit$T)

    y.star <- sapply(1:ncol(mf$B), function(i) 
    {
	grp <- (!!fit$Z[ , i])[groups]
	mu <- newdata[ , grp] %*% fit$B[grp, i]
	rnorm(n, mu, sqrt(tau))
    })

    y.star <- matrix(y.star, nrow=n)

    res <- list(
	mean=apply(y.star, 1, mean),
	quantiles=apply(y.star, 1, quantile, probs=quantiles)
    )

    if (return_samples) {
	res <- c(res, list(samples=y.star))
    }

    return(res)
}


pred <- spsl.predict(mf, d$groups, d$X, return_samples=T)
# coverage
mean(
    (pred$quantiles[1, ] <= d$y) & (pred$quantiles[2, ] >= d$y)
)

pred <- spsl.predict(mf, d$groups, matrix(d$X[1, ], nrow=1), return_samples=T)
plot(density(pred$samples))
abline(v=d$y[1])


# ----------------------------------------
# Comparison
# ----------------------------------------
gpred <- gsvb.predict(f, d$groups, matrix(d$X[1, ], nrow=1), 1e5, return_samples=T)
spred <- spsl.predict(mf, d$groups, matrix(d$X[1, ], nrow=1), return_samples=T)

pdf("~/p1/notes/gsvb/figures/post_pred.pdf", width=6, height=4)
    par(mar=c(2,4,1,0))
    plot(density(gpred$samples[1, ]), lwd=2, col=2, main="")
    lines(density(spred$samples[1, ]), lwd=2, col=3, lty=2)
    grid()
    legend("topright", legend=c("Variational PP", "PP"), lwd=2, col=c(2,3), lty=c(1,2))
dev.off()


gpred <- gsvb.predict(f,  d$groups,  d$X, 1e3, return_samples=T)
spred <- spsl.predict(mf, d$groups, d$X, return_samples=T)
mean(
    (gpred$quantiles[1, ] <= d$y) & (gpred$quantiles[2, ] >= d$y)
)
mean(
    (spred$quantiles[1, ] <= d$y) & (spred$quantiles[2, ] >= d$y)
)


# ----------------------------------------
# Testing dgp with test data
# ----------------------------------------
source("00-functions.R")
n <- 200
p <- 1000 
g <- 5
s <- 3

d <- dgp_block(200, 1000, 5, 3, list(corr=0.6, block_size=50), seed=91)
d$test


# ----------------------------------------
# testing method_post_pred
# ----------------------------------------
d <- dgp_block(200, 1000, 5, 3, list(corr=0.6, block_size=50), seed=91)

m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3, diag_covariance=TRUE)
f <- m_gsvb(d, m_par=m_par)

m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3, diag_covariance=FALSE)
f <- m_gsvb(d, m_par=m_par)

f <- gsvb::gsvb.fit(d$y, d$X, d$groups, diag_covariance=F)
(coverage <- method_post_pred(d, f, "gsvb"))


# ----------------------------------------
# testing method coverage
# ----------------------------------------
d <- dgp_block(200, 1000, 5, 3, list(corr=0.6, block_size=50), seed=91)

m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3, diag_covariance=TRUE)
f <- m_gsvb(d, m_par=m_par)

m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3, diag_covariance=FALSE)
f <- m_gsvb(d, m_par=m_par)

f <- gsvb::gsvb.fit(d$y, d$X, d$groups, diag_covariance=F)
method_coverage(d, f, "gsvb")

s <- spsl::spsl.fit(d$y, d$X, d$groups)
method_coverage(d, s, "spsl")

