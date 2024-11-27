# Misc code for debugging
source("./00-functions.R")

# ----------------------------------------
# SpSL MCMC
# ----------------------------------------
d <- dgp_wishart(200, 1000, 10, 10, 1.5, list(model="gaussian", corr=3, weight=0.9), seed=2)

f <- spsl::spsl.fit(d$y, d$X, d$groups, mcmc_samples=5e4, burnin=3e4, a_0 = 1, b_0=200)

plot(f$beta_hat)
points(d$b, pch=20)
plot(f$Z[35, ])
plot(f$Z[40, ])
plot(f$Z[2, ])
points(c(0, d$b), pch=20)
plot(c(0, d$b), pch=20)
plot(f$g)

f0 <- gsvb::gsvb.fit(d$y, d$X, d$groups, diag_covariance=FALSE, niter=500)
plot(f0$beta_hat)
points(c(0, d$b), pch=20)
plot(f0$g)

# length(unique(d$groups))

w <- rbeta(1, 190, 200)
z <- f$Z[ , 5000][-1]
z.old <- z[d$groups]
z <- z.old
b <- f$B[ , 5000][-1]


vals <- c()
for (i in 1:100)
{
    G <- which(d$groups == i)

    z[G] <- 0
    lp0 = sum(dnorm(d$y, d$X[ , which(!!z)] %*% d$b[which(!!z)], 1, log=T)) + dbeta(1-w, 1, 200, T)

    z[G] <- 1
    lp1 = sum(dnorm(d$y, d$X[ , which(!!z)] %*% d$b[which(!!z)], 1, log=T)) + dbeta(1-w, 1, 200, T)
    
    z[G] <- z.old[G]
    vals <- c(vals, print(1/(1+exp(lp0 - lp1))))
}
vals[d$active]


dnorm(y, X %*% f$b, 1, log=T)


matplot(t(f$B[which(!!d$b)[1:10], ]), type="l")
matplot(t(f$Z[which(!!d$b)[1:10], ]), type="l")
matplot(t(f$B[which(!!d$b)[1:10], ]), type="l")
matplot(t(f$B[which(!!d$b)[11:20], ]), type="l")
matplot(t(f$B[which(!!d$b)[12:20], ]), type="l")
plot(f$B[which(!!d$b)[1], ], type="l")
plot(f$B[which(!!d$b)[2], ], type="l")
plot(f$B[which(!!d$b)[1], ], type="l")
plot(f$B[which(!!d$b)[1], ], type="l")
plot(f$B[which(!!d$b)[1], ], type="l")
plot(f$B[which(!!d$b)[21], ], type="l")
plot(f$B[which(!!d$b)[22], ], type="l")
plot(f$B[which(!!d$b)[23], ], type="l")
plot(f$B[which(!!d$b)[24], ], type="l")


newdata <- d$test$X
newdata <- cbind(1, newdata)

# mu <- sapply(1:4500, function(i) 
# {
#     grp <- (!!f$Z[ , i])[f$parameters$groups]
#     if (any(grp)) {
# 	if (sum(grp) == 1) {
# 	    xb <- newdata[ , grp] * f$B[grp, i] 
# 	} else {
# 	    xb <- newdata[ , grp] %*% f$B[grp, i]
# 	}
#     } else {
# 	xb <- rep(0, nrow(newdata))
#     }
#     return(xb)
# })


plot(d$test$y)
points(pred$mean, pch=20)
points(pred$quantiles[1, ])
points(pred$quantiles[2, ])

pred <- spsl::spsl.predict(f, newdata=d$test$X)

# ----------------------------------------
# Bimom settings
# ----------------------------------------
d <- dgp_diag(1000, 5000, 5, 3, 1.5, list(model="binomial", corr=0))


d <- dgp_wishart(350, 1000, 5, 3, 0.8, list(model="binomial", dof=3, weight=0.9))
m_par <- list(family="binomial-jensens", lambda=1, a0=1, b0=200, diag_covariance=FALSE, intercept=TRUE)
m_par <- list(family="binomial-jaakkola", lambda=1, a0=1, b0=200, diag_covariance=TRUE, intercept=TRUE)
m_par <- list(family="binomial-jaakkola", lambda=1, a0=1, b0=200, diag_covariance=FALSE, intercept=TRUE)



d <- dgp_wishart(400, 1000, 5, 3, 0.8, list(model="binomial", dof=3, weight=0.9))
f <- spsl::spsl.fit(d$y, d$X, d$groups, family="binomial", 
		    mcmc_samples=1e4, burnin=5e3)


f <- m_gsvb(d, m_par)

d <- dgp_wishart(500, 5000, 10, 10, 0.8, list(model="binomial", dof=3, weight=0.9))
m_par <- list(family="binomial", l0=100, l1=1, a0=1, b0=5000/5)
m_par <- list(family="binomial", l0=20, l1=1, a0=1, b0=5000/5)
f <- m_ssgl(d, m_par)

d <- dgp_wishart(350, 1000, 5, 5, 1.5, list(model="binomial", dof=3, weight=0.9))


m_par <- list(family="binomial", lambda=1, a0=1, b0=200, mcmc_samples=1e5, burnin=5e4, intercept=TRUE)


f <- m_spsl(d, m_par)

d <- dgp_diag(350, 1000, 5, 5, 1.5, list(model="binomial", corr=0))
d <- dgp_wishart(350, 1000, 5, 5, 1.5, list(model="binomial", dof=3, weight=0.9))
f <- spsl::spsl.fit(d$y, d$X, d$groups, family="binomial", 
		    mcmc_samples=1e4, burnin=5e3)
plot(f$beta)
matplot(t(f$B[ c(F, !!d$b), ]), type="l")
points(d$b, pch=20)

d <- dgp_diag(500, 5000, 5, 5, 0.8, list(model="binomial", corr=0))




d <- dgp_diag(500, 2500, 5, 3, 1.0, list(model="binomial", corr=0))
d <- dgp_diag(500, 5000, 5, 3, 1.0, list(model="binomial", corr=0))
d <- dgp_wishart(500, 2500, 5, 3, 1.0, list(model="binomial", dof=3, weight=0.9))

f <- gsvb::gsvb.fit(d$y, d$X, d$groups, family="binomial-jaakkola")

plot(f$beta)
points(d$b, pch=20)


# ----------------------------------------
# Pois settings
# ----------------------------------------
d <- dgp_wishart(500, 1000, 5, 10, 0.45, list(model="poisson", dof=3, weight=0.9), seed=5)

m_par <- list(family="poisson", lambda=1, a0=1, b0=200, 
	      diag_covariance=TRUE, intercept=FALSE)

f <- m_gsvb(d, m_par)
f
f <- gsvb::gsvb.fit(d$y, d$X, d$groups, family="poisson", intercept=FALSE, 
		    diag_covariance=TRUE, niter=250)
dev.off()
plot(f$b)
points(d$b, pch=20)


f <- gsvb::gsvb.fit(d$y, d$X, d$groups, family="poisson", intercept=FALSE, diag_covariance=FALSE,
		    s=rep(0.25, 1000), niter=35)

plot(f$beta)
points(d$b, pch=20)

plot(d$b, pch=20)
points(f$beta)
plot(f$s^2)
points(d$b, pch=20)


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
d <- dgp_block(200, 1000, 5, 3,list(corr=0.6, block_size=50), seed=91)
dgp_block()

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


# ----------------------------------------
# testing mcmc
# ----------------------------------------
d <- dgp_diag(200, 1000, 10, 10, 1.5, list(model="gaussian", corr=0), seed=1)
d <- dgp_diag(200, 1000, 5, 5, 1.5, list(model="gaussian", corr=0), seed=1)

d <- dgp_diag(300, 500, 5, 5, 1.5, list(model="gaussian", corr=0), seed=1)
f0 <- spsl::spsl.fit(d$y, d$X, d$groups, family="gaussian", intercept=T, 
		    mcmc_samples=20e3, burnin = 5e3, thin=1)

matplot(t(f0$B[1:10, ]), type="l")
plot(f0$g)
plot(f0$beta_hat)
points(d$b, pch=20)
method_coverage(d, f0, "spsl", 0.95)


mean(
    d$b[!!d$b] >=
    apply(f0$B[1 + which(d$b != 0), ], 1, quantile, prob=c(0.025, 0.975))[1, ] &
    d$b[!!d$b] <=
    apply(f0$B[1 + which(d$b != 0), ], 1, quantile, prob=c(0.025, 0.975))[2, ]
)

gsvb::gsvb.credible_intervals

plot(d$b[!!d$b], pch=20)
points(spsl::spsl.credible_intervals(f0, 0.99)[1 + which(d$b != 0) , 1])
points(spsl::spsl.credible_intervals(f0, 0.99)[1 + which(d$b != 0) , 2])

# f1 <- spsl::spsl.fit(d$y, d$X, d$groups, family="gaussian", intercept=T, 
# 		    mcmc_samples=20e3, burnin = 10e3)
# f2 <- spsl::spsl.fit(d$y, d$X, d$groups, family="gaussian", intercept=T, 
# 		    mcmc_samples=20e3, burnin = 10e3)
# f3 <- spsl::spsl.fit(d$y, d$X, d$groups, family="gaussian", intercept=T, 
# 		    mcmc_samples=20e3, burnin = 10e3)

spsl::spsl.credible_intervals
method_coverage(d, f0, "spsl", 0.95)
plot(d$b[!!d$b], pch=20)
points(spsl::spsl.credible_intervals(f0, 0.99)[1 + which(d$b != 0) , 1])
points(spsl::spsl.credible_intervals(f0, 0.99)[1 + which(d$b != 0) , 2])

f1 = gsvb::gsvb.fit(d$y, d$X, d$groups, diag_covariance=FALSE)
f1 = gsvb::gsvb.fit(d$y, d$X, d$groups, diag_covariance=TRUE)
method_coverage(d, f1, "gsvb", prob=0.95)
gsvb::gsvb.credible_intervals


matplot(t(f0$B[f0$parameters$groups == 44, ]), type="l")
d$b[45 * 5 + 1:5]
points(d$b[d$groups == 43])

l = coda::mcmc.list(
	coda::mcmc(t(f0$B)),
	coda::mcmc(t(f1$B)),
	coda::mcmc(t(f2$B)),
	coda::mcmc(t(f3$B))
    )

perf = coda::gelman.diag(l)
plot(perf$psrf[ , 1])
perf$mpsrf
plot(f0$g)
plot(f1$g)
plot(f2$g)
plot(f3$g)





plot(coda::effectiveSize(l))

plot(f0$g, ylim=c(0, 1))
points(f1$g)
points(f2$g)
points(f3$g)


install.packages("coda")

coda::gelman.diag
coda::mcmc.list
mean(d$b[!!d$b])

layout(1)
matplot(t(f$B[c(FALSE, !!d$b), ]), type="l")



d <- dgp_diag(200, 1000, 10, 10, 1.5, list(model="gaussian", corr=0), seed=1)
d <- dgp_diag(200, 1000, 5, 5, 1.5, list(model="gaussian", corr=0), seed=1)
m_par <- list(family="gaussian", lambda=1, a0=1, b0=1000/5 + 1, a_t=1e-3, 
	      b_t=1e-3, mcmc_samples=2e3, burnin=5e2, intercept=TRUE)
f <- m_spsl(d, m_par)

d <- dgp_diag(400, 1000, 5, 3, 0.45, list(model="poisson", corr=0), seed=79)
f = spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", intercept=FALSE,
    mcmc_samples=100000, burnin=50000)

f = spsl::spsl.fit(d$y, d$X, d$groups)
f$B * 


f$B * f$Z[f$parameters$groups, ]



plot(f$beta_hat)
points(d$b, pch=20)

matplot(t(f$B[1:10,]), type="l")

setwd("./r/simulations")
source("./00-functions.R")
d <- dgp_diag(200, 1000, 10, 10, 1.5, list(model="gaussian", corr=0), seed=1)
p = 1000 
g = 200
m_par = list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
    diag_covariance=TRUE, intercept=TRUE, ordering=0, init_method="lasso")
    
m_par = list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
    diag_covariance=TRUE, intercept=TRUE, ordering=0, init_method="ridge")

m_par = list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
    diag_covariance=TRUE, intercept=TRUE, ordering=0)

ff = m_gsvb(d, m_par)

f = gsvb::gsvb.fit(d$y, d$X, d$groups,
diag_covariance = TRUE, intercept = TRUE, ordering=0)

plot(f$mu)


p = 1000 
g = 200
# m_par <- list(family="gaussian", lambda=1, a0=1, b0=1000/5 + 1, a_t=1e-3, 
# 	      b_t=1e-3, mcmc_samples=1e3, burnin=0, intercept=TRUE)
setwd("./r/simulations")
source("./00-functions.R")



# =============================================================================
#
#                               GAUSSIAN
#
# =============================================================================
library(ParBayesianOptimization)

# Define the objective function for Bayesian Optimization
objective_function <- function(k1, k2) {
    seed = sample(1:100, 1)
  d <- dgp_diag(200, 1000, 5, 5, bmax=1.5, list(model="gaussian", corr=0.6), seed=seed)
  
  # Fit the model with the given kernel parameters
  f0 <- spsl::spsl.fit(d$y, d$X, d$groups, family="gaussian", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  f1 <- spsl::spsl.fit(d$y, d$X, d$groups, family="gaussian", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  f2 <- spsl::spsl.fit(d$y, d$X, d$groups, family="gaussian", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  f3 <- spsl::spsl.fit(d$y, d$X, d$groups, family="gaussian", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  
  # Create MCMC list
  l <- coda::mcmc.list(
    coda::mcmc(t(f0$B)),
    coda::mcmc(t(f1$B)),
    coda::mcmc(t(f2$B)),
    coda::mcmc(t(f3$B))
  )
  
  # Calculate Gelman-Rubin diagnostic
  tryCatch({
    gg <- coda::gelman.diag(l, multivariate=TRUE)
    
    # Return the multivariate potential scale reduction factor (mpsrf) as the objective
    return(list(Score = -gg$mpsrf))  # Negative because we want to minimize mpsrf
  }, error = function(e) {
    cat("Error in calculating Gelman-Rubin diagnostic: ", e$message, "\n")
    return(list(Score = -500))
  })
}

# Run Bayesian Optimization
opt_results <- bayesOpt(
  FUN = objective_function,
  bounds = list(k1 = c(0.05, 0.2), k2 = c(10, 25)),
  initPoints = 5,
  iters.n = 15,
  acq = "ei"
)

# Print the optimal parameters
opt_results$scoreSummary
ParBayesianOptimization::getBestPars(opt_results)

objective_function(0.15, 19.5)



kernel_param_1=0.05, kernel_param_2=10
# these 
# these are good for poisson
kernel_param_1=0.025, kernel_param_2=20





# =============================================================================
#
#                               POISSON
#
# =============================================================================
objective_function <- function(k1, k2) {
  seed = sample(1:100, 1)
  d <- dgp_diag(400, 1000, 5, 2, bmax=0.45, list(model="poisson", corr=0.6), seed=seed)
  
  # Fit the model with the given kernel parameters
  f0 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  f1 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  f2 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  f3 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  
  # Create MCMC list
  l <- coda::mcmc.list(
    coda::mcmc(t(f0$B)),
    coda::mcmc(t(f1$B)),
    coda::mcmc(t(f2$B)),
    coda::mcmc(t(f3$B))
  )
  
  # Calculate Gelman-Rubin diagnostic
  tryCatch({
    gg <- coda::gelman.diag(l, multivariate=TRUE)
    
    # Return the multivariate potential scale reduction factor (mpsrf) as the objective
    return(list(Score = -gg$mpsrf))  # Negative because we want to minimize mpsrf
  }, error = function(e) {
    cat("Error in calculating Gelman-Rubin diagnostic: ", e$message, "\n")
    return(list(Score = -500))
  })
}

opt_results_pois <- bayesOpt(
  FUN = objective_function,
  bounds = list(k1 = c(0.01, 0.1), k2 = c(15, 25)),
  initPoints = 5,
  iters.n = 20,
  acq = "ei"
)

# Print the optimal parameters
opt_results_pois$scoreSummary
ParBayesianOptimization::getBestPars(opt_results_pois)




# =============================================================================
#
#                               MISC
#
# =============================================================================

f0 = spsl::spsl.fit(d$y, d$X, d$groups, family="binomial", mcmc_samples = 1.5e4, burnin = 5e3, kernel_param_1=0.2, kernel_param_2=9)
f1 = spsl::spsl.fit(d$y, d$X, d$groups, family="binomial", mcmc_samples = 1.5e4, burnin = 5e3, kernel_param_1=0.2, kernel_param_2=9)
f2 = spsl::spsl.fit(d$y, d$X, d$groups, family="binomial", mcmc_samples = 1.5e4, burnin = 5e3, kernel_param_1=0.2, kernel_param_2=9)
f3 = spsl::spsl.fit(d$y, d$X, d$groups, family="binomial", mcmc_samples = 1.5e4, burnin = 5e3, kernel_param_1=0.2, kernel_param_2=9)

l = coda::mcmc.list(
coda::mcmc(t(f0$B # * f0$Z[f0$parameters$groups, ]	
)),
coda::mcmc(t(f1$B # * f1$Z[f1$parameters$groups, ]	
)),
coda::mcmc(t(f2$B # * f2$Z[f2$parameters$groups, ]	
)),
coda::mcmc(t(f3$B # * f3$Z[f3$parameters$groups, ]	
)))

gg = coda::gelman.diag(l, multivariate = TRUE)
quantile(gg$psrf[,2], 0.95, na.rm=TRUE)
gg$mpsrf
# hist(gg$psrf[1, ])

f0$g
d$active_groups
which(f0$parameters$groups == (d$active_groups + 1)[2])

plot(f0$g)
sum(f0$g > 0.5)
plot(f0$beta_hat)

plot(f0$beta_hat)
points(f1$beta_hat, col="red")
points(f2$beta_hat, col="green")
points(f3$beta_hat, col="blue")

f0$g
plot(f0$B[10, ], type="l")
plot(f0$B[1, ], type="l")
plot(f0$B[2, ], type="l")
plot(f0$B[3, ], type="l")
plot(f0$B[4, ], type="l")
plot(f0$B[2, ], type="l")

d$active_groups

colSums(fit$Z) 
# gg$mpsrf
# hist(gg$psrf[ , 1])
# gg$mpsrf

matplot( t(f3$B * f3$Z[f3$parameters$groups, ] ), type="l")
plot(f3$B[78, ], type="l")

library(coda)
gg = my.gelman.diag(l)

chol(gg$W)


min(eigen(gg$W, only.values = TRUE)$values)

matrixNormal::is.positive.definite(gg$W, tol = 1e-18)


