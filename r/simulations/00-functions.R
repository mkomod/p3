library(mvtnorm)


# ----------------------------------------
# Data generating processes
# ----------------------------------------
.dgp_base <- function(n, p, gsize, s, b, seed)
{
    set.seed(seed)

    if (sum(gsize) == p) {
	groups <- rep(1:length(gsize), each=gsize)
    } else {
	if (p %% gsize != 0) stop("number of groups must be a factor of p")
	groups <- rep(1:(p/gsize), each=gsize)
    }
    
    active_groups <- sample(1:length(unique(groups)), s)

    if (is.null(b)) {
	b <- rep(0, p)
	for (a in active_groups) {
	    gk <- which(groups == a)
	    m <- length(gk)
	    b[gk] <- sample(c(1, -1), m, replace=T) * runif(m, min=1, max=3)
	}
    }
    
    return(list(groups=groups, active_groups=active_groups, b=b))
}


.dgp_return <- function(n, n.test, y, X) 
{
    if (n.test == 0) {
	return(list(y=y, X=X))
    } else {
	N <- n + n.test
	return(list(y=y[1:n], X=X[1:n, ], 
		    test=list(y=y[(n+1):N], X=X[(n+1):N, ])))
    }
}


dgp_diag <- function(n, p, gsize, s, pars, b=NULL, seed=1, sig=1, n.test=100)
{
    N <- ifelse(n.test > 0, n + n.test, n)
    res <- .dgp_base(N, p, gsize, s, b, seed)

    # unpack pars
    corr <- pars[[1]]

    if (corr == 0) {
	X <- matrix(rnorm(N * p), nrow=N)
    } else {
	S <- outer(1:p, 1:p, function(i, j) corr^abs(i - j))
	X <- mvtnorm::rmvnorm(N, rep(0, p), sigma=S)
    }
    
    j <- which(res$groups %in% res$active_groups)
    y <- X[ , j] %*% res$b[j] + rnorm(N, sd=sig)
    
    dat <- .dgp_return(n, n.test, y, X)
    return(c(res, dat, list(seed=seed, sig=sig, corr=corr)))
}


dgp_block <- function(n, p, gsize, s, pars, b=NULL, seed=1, sig=1, n.test=100)
{
    N <- ifelse(n.test > 0, n + n.test, n)
    res <- .dgp_base(N, p, gsize, s, b, seed)

    # unpack pars
    corr <- pars[[1]]
    block_size <- pars[[2]]

    if (corr == 0) {
	X <- matrix(rnorm(N * p), nrow=N)
    } else {
	X <- matrix(nrow=N, ncol=0)
	for (block in 1:(p / block_size)) {
	    S <- matrix(corr, nrow=block_size, ncol=block_size)
	    diag(S) <- 1
	    X <- cbind(X, mvtnorm::rmvnorm(N, mean=rep(0, block_size), S))
	}
    }

    j <- which(res$groups %in% res$active_groups)
    y <- X[ , j] %*% res$b[j] + rnorm(N, sd=sig)

    dat <- .dgp_return(n, n.test, y, X)
    return(c(res, dat, list(seed=seed, sig=sig, corr=corr)))
}


dgp_wishart <- function(n, p, gsize, s, pars=list(dof=3, weight=0.9), 
    b=NULL, seed=1, sig=1, n.test=100)
{
    N <- ifelse(n.test > 0, n + n.test, n)
    res <- .dgp_base(N, p, gsize, s, b, seed)

    # unpack pars
    dof <- pars[[1]]
    weight <- pars[[2]]

    S <- (1-weight)*solve(rWishart(1, p+dof, diag(rep(1, p)))[, , 1])
    for (block in 1:(p / gsize)) {
	# generate cov for group
	S_g <- solve(rWishart(1, gsize+dof, diag(rep(1, gsize)))[, , 1])

	G <- (1+(block-1)*gsize):(gsize*block)
	S[G, G] <- S[G, G] + weight*S_g
    }
    X <- mvtnorm::rmvnorm(N, rep(0, p), sigma=S)

    j <- which(res$groups %in% res$active_groups)
    y <- X[ , j] %*% res$b[j] + rnorm(N, sd=sig)

    dat <- .dgp_return(n, n.test, y, X)
    return(c(res, dat, list(seed=seed, sig=sig, dof=dof, weight=weight)))
}


# ----------------------------------------
# Methods (parallelized)
# ----------------------------------------
m_run <- function(method, method_par, setting_par, CORES) 
{
    # create a new env for cluster
    e <- list2env(setting_par)
    e$method_par <- method_par
    e$.dgp_base <- .dgp_base
    e$.dgp_return <- .dgp_return
    e$method <- method
    e$method_summary <- method_summary
    e$method_post_pred <- method_post_pred
    
    # init cluster
    cl <- parallel::makeCluster(getOption("cl.cores", CORES), outfile="")
    parallel::clusterExport(cl, ls(e, all.names=TRUE), envir=e)
    on.exit(parallel::stopCluster(cl))
    
    # run the method
    res <- parallel::parSapply(cl, 1:runs, function(run) {
	cat(run, "\t")
	d <- dgp(n, p, g, s, pars, seed=run)
	method(d, method_par)
    })

    return(t(res))
}


m_gsvb <- function(d, m_par=list(lambda=0.5, a0=1, b0=100, a_t=1e-3, b_t=1e-3,
	diag_covariance=TRUE))
{
    tryCatch({
	fit.time <- system.time({
	    fit <- gsvb::gsvb.fit(d$y, d$X, d$groups, intercept=TRUE, 
		diag_covariance=m_par$diag_covariance, lambda=m_par$lambda,
		a0=m_par$a0, b0=m_par$b0, tau_a0=m_par$a_t, tau_b0=m_par$b_t,
		niter=500, track_elbo=FALSE)
	})

	active_groups <- rep(0, length(unique(d$groups)))
	active_groups[d$active_groups] <- 1
	res <- method_summary(d$b, active_groups, fit$beta_hat[-1], fit$g[-1], 0.5)
	
	if (!is.null(d$test)) {
	    coverage <- method_post_pred(d, fit, method="gsvb", 
		quantiles=c(0.025, 0.975), return_samples=FALSE)

	    return(c(unlist(res), unlist(fit.time[3]), unlist(coverage)))
	}

	return(c(unlist(res), unlist(fit.time[3])))

    }, error=function(e) 
    {
	print(e)
	cat("error in run: ", d$seed)
	return(rep(NA, 214))
    })
}


m_spsl <- function(d, m_par=list(family="linear", lambda=0.5, a0=1, b0=100, 
    a_t=1e-3, b_t=1e-3, mcmc_samples=10e3))
{
    fit.time <- system.time({
	fit <- spsl::spsl.group_sparse(d$y, d$X, d$groups, family=m_par$family,
	    lambda=m_par$lambda, a_0=m_par$a0, b_0=m_par$b0, a_t=m_par$a_t,
	    b_t=m_par$b_t, mcmc_sample=m_par$mcmc_samples)
    })

    active_groups <- rep(0, length(unique(d$groups)))
    active_groups[d$active_groups] <- 1
    res <- method_summary(d$b, active_groups, fit$beta_hat[-1], fit$g[-1], 0.5)

    if (!is.null(d$test)) {
	coverage <- method_post_pred(d, fit, method="spsl", 
	    quantiles=c(0.025, 0.975), return_samples=FALSE)

	return(c(unlist(res), unlist(fit.time[3]), unlist(coverage)))
    }

    return(c(unlist(res), unlist(fit.time[3])))
}


m_ssgl <- function(d, m_par=list(l0=20, l1=1, a0=1, b0=100)) 
{
    fit.time <- system.time({
	fit <- sparseGAM::SSGL(d$y, d$X, d$X, d$groups, family="gaussian",
	    lambda0=m_par$l0, lambda1=m_par$l1, a=m_par$a0, b=m_par$b0)
    })

    active_groups <- rep(0, length(unique(d$groups)))
    active_groups[d$active_groups] <- 1
    res <- method_summary(d$b, active_groups, fit$beta, fit$classifications, 0.5)
    return(c(unlist(res), unlist(fit.time[3])))
}


# ----------------------------------------
# Method evaluation
# ----------------------------------------
method_summary <- function(beta_true, active_groups, beta_hat, 
    inclusion_prob, threshold) 
{
    g <- beta_true != 0
    tpr <- numeric(0)
    fpr <- numeric(0)

    for (thresh in seq(0, 1, by=0.01)) {
	tab <- table(inclusion_prob >= thresh, active_groups)
	if (nrow(tab) == 2) {
	    TP <- tab[2,2]; FP <- tab[2,1]
	    FN <- tab[1,2]; TN <- tab[1,1]
	} else if(all(inclusion_prob >= thresh)) {
	    TP <- sum(g);   FP <- sum(1 - g)
	    FN <- 0; 	    TN <- 0
	} else {
	    TP <- 0;        FP <- 0
	    FN <- sum(g);   TN <- sum(1 - g)
	}
	tpr <- c(tpr, TP / (TP + FN))
	fpr <- c(fpr, FP / (TN + FP))
    }

    # compute the AUC using trapizium int
    auc <- -sum((tpr[1:100] + tpr[2:101]) / 2 * diff(fpr))


    tab <- table(inclusion_prob > threshold, active_groups)
    if (nrow(tab) == 2) {
	TP <- tab[2,2]; FP <- tab[2,1]
	FN <- tab[1,2]; TN <- tab[1,1]
    } else if(all(inclusion_prob > thresh)) {
	TP <- sum(g);   FP <- sum(1 - g)
	FN <- 0;        TN <- 0
    } else {
	TP <- 0;        FP <- 0
	FN <- sum(g);   TN <- sum(1 - g)
    }

    return(list(
	acc = (TN + TP) / (TN + TP + FN + FP), 
	err = (FP + FN) / (TN + TP + FN + FP), 
	tpr = TP / (TP + FN), 
	tnr = TN / (TN + FP), 
	fpr = FP / (TN + FP),
	ppv = ifelse(TP + FP == 0, NA, TP / (TP + FP)),
	fdr = ifelse(TP + FP == 0, NA, FP / (TP + FP)),
	auc=auc,
	FPR=fpr,
	TPR=tpr,
	l1 = sum(abs(beta_true - beta_hat)),
	l2 = sqrt(sum((beta_true - beta_hat)^2))
    ))
}


method_post_pred <- function(d, fit, method, quantiles=c(0.025, 0.975), 
	return_samples=FALSE, mcn=1e4) 
{
    if (method == "spsl") {
	pp <- spsl::spsl.predict(fit, d$test$X, quantiles=quantiles,
	    return_samples=return_samples)
    } 
    if (method == "gsvb") {
	pp <- gsvb::gsvb.predict(fit, d$test$X, mcn=mcn, 
	    quantiles=quantiles, return_samples=return_samples)
    }
    coverage <- mean((pp$quantiles[1, ] <= d$test$y) & (d$y <= pp$quantiles[2, ]))

    return(list(coverage=coverage))
}


method_coverage <- function()
{
    stop("not implemented")
}


# Compute the coverage
# method_coverage <- function(fit, d, a=0.05, threshold=0.5, conditional=FALSE, 
# 	mcmc=FALSE)
# {
#     if (mcmc) {
# 	ci <- mcmc.credible_interval(fit, a, conditional=conditional)
#     } else {
# 	ci <- svb.credible_interval(fit, a, conditional=conditional)
#     }
#     l <- ci[ , 1]
#     u <- ci[ , 2]
#     dirac <- ci[ , 3]

#     coverage <- (((l <= d$beta) & (d$beta <= u)) | ((d$beta == 0) & (dirac == T)))
#     bs <- d$beta != 0

#     return(list(
# 	coverage.n0 = mean(coverage[bs]),
# 	length.n0 = mean(u[bs] - l[bs]),
# 	coverage.0 = mean(coverage[!bs]),
# 	length.0 =  mean(u[!bs] - l[!bs])
#     ))
# }


# ----------------------------------------
# Misc
# ----------------------------------------
# import an environment variable ex: cores <- d("CORES", 8)
read.env <- function(evar, default) 
{
    if(Sys.getenv(evar) != "") as(Sys.getenv(evar), class(default)) else default
}


# # compute the credible intervals
# mcmc.credible_interval <- function(fit, a=0.05, conditional=FALSE, burnin=1e3)
# {
#     p <- ncol(fit$b)

#     credible.interval <- sapply(1:length(fit$g), function(i) {
# 	g <- fit$g[i]	

# 	if (conditional) {
# 	    m <- fit$b[i, burnin:p]
# 	    z <- !!fit$z[i, burnin:p]
# 	    m <- m[z]
# 	    f <- density(m)
# 	    cdf <- ecdf(m)
# 	    return(quantile(cdf, c(a/2, 1 - a/2)))
# 	}

# 	if (!conditional) {
# 	    if (g > 1 - a) {
# 		# Slab contains 1 - a of mass
# 		a.g <- 1 - (1-a)/g
# 		m <- fit$b[i, burnin:p]
# 		z <- !!fit$z[i, burnin:p]
# 		m <- m[z]
# 		f <- density(m)
# 		cdf <- ecdf(m)
# 		interval <- quantile(cdf, c(a.g/2, 1 - a.g/2))
# 		contains.dirac <- FALSE
		
# 		if (interval[1] <= 0 && interval[2] >= 0) {
# 		    # if interval contains Dirac mass it needs to be smaller
# 		    interval <- quantile(cdf, c(a.g/2+(1-g)/2, 1-a.g/2-(1-g)/2))
# 		    contins.dirac <- TRUE
# 		}

# 		return(c(lower=interval[1], upper=interval[2], 
# 			 contains.dirac=contains.dirac))
# 	    } else if (g < a) {
# 		# Dirac contains 1 - a of mass
# 		return(c(lower=0, upper=0, contains.dirac=T))
# 	    } else {
# 		# will always contain the Dirac
# 		m <- fit$b[i, burnin:p]
# 		z <- !!fit$z[i, burnin:p]
# 		m <- m[z]
# 		f <- density(m)
# 		cdf <- ecdf(m)
# 		interval <- quantile(cdf, c(a/2 +(1-g)/2, 1 - a/2 - (1-g)/2))

# 		return(c(lower=interval[1], upper=interval[2], 
# 			 contains.dirac=TRUE))
# 	    }
# 	}
#     })

#     return(t(credible.interval))
# }
