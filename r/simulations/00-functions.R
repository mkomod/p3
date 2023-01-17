library(mvtnorm)


# ----------------------------------------
# Data generating processes
# ----------------------------------------
.dgp_base <- function(n, p, gsize, s, bmax, b, seed)
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
	    b[gk] <- sample(c(1, -1), m, replace=T) * runif(m, min=0, max=bmax)
	}
    }
    
    return(list(groups=groups, active_groups=active_groups, b=b))
}


.dgp_return <- function(model, X, Xb, n, n.test, sig=1)
{
    if (model == "gaussian") {	# LINEAR REG
	y <- Xb + rnorm(n + n.test, 0, sd=sig)
    }
    if (model == "binomial") {      # LOGISTIC REG
	probs <- 1/(1+exp(-Xb))
	y <- rbinom(n + n.test, 1, probs)
    }
    if (model == "poisson") {	# POISSON REG
	y <- rpois(n + n.test, lambda=exp(Xb))
    }

    if (n.test == 0) {
	return(list(y=y, X=X))
    } else {
	N <- n + n.test
	return(list(y=y[1:n], X=X[1:n, ], 
		    test=list(y=y[(n+1):N], X=X[(n+1):N, ])))
    }
}


dgp_diag <- function(n, p, gsize, s, bmax, pars, b=NULL, seed=1, sig=1, n.test=100)
{
    N <- ifelse(n.test > 0, n + n.test, n)
    res <- .dgp_base(N, p, gsize, s, bmax, b, seed)

    # unpack pars
    model <- pars[[1]]
    corr <- pars[[2]]

    if (corr == 0) {
	X <- matrix(rnorm(N * p), nrow=N)
    } else {
	S <- outer(1:p, 1:p, function(i, j) corr^abs(i - j))
	X <- mvtnorm::rmvnorm(N, rep(0, p), sigma=S)
    }
    
    j <- which(res$groups %in% res$active_groups)
    Xb <- X[ , j] %*% res$b[j]
    
    dat <- .dgp_return(model, X, Xb, n, n.test, sig)
    return(c(res, dat, list(seed=seed, sig=sig, corr=corr)))
}


dgp_block <- function(n, p, gsize, s, bmax, pars, b=NULL, seed=1, sig=1, n.test=100)
{
    N <- ifelse(n.test > 0, n + n.test, n)
    res <- .dgp_base(N, p, gsize, s, bmax, b, seed)

    # unpack pars
    model <- pars[[1]]
    corr <- pars[[2]]
    block_size <- pars[[3]]

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
    Xb <- X[ , j] %*% res$b[j]

    dat <- .dgp_return(model, X, Xb, n, n.test, sig)
    return(c(res, dat, list(seed=seed, sig=sig, corr=corr)))
}


dgp_wishart <- function(n, p, gsize, s, bmax, pars=list(dof=3, weight=0.9), 
    b=NULL, seed=1, sig=1, n.test=100)
{
    N <- ifelse(n.test > 0, n + n.test, n)
    res <- .dgp_base(N, p, gsize, s, bmax, b, seed)

    # unpack pars
    model <- pars[[1]]
    dof <- pars[[2]]
    weight <- pars[[3]]

    S <- (1-weight)*solve(rWishart(1, p+dof, diag(rep(1, p)))[, , 1])
    for (block in 1:(p / gsize)) {
	# generate cov for group
	S_g <- solve(rWishart(1, gsize+dof, diag(rep(1, gsize)))[, , 1])

	G <- (1+(block-1)*gsize):(gsize*block)
	S[G, G] <- S[G, G] + weight*S_g
    }
    X <- mvtnorm::rmvnorm(N, rep(0, p), sigma=S)

    j <- which(res$groups %in% res$active_groups)
    Xb <- X[ , j] %*% res$b[j]

    dat <- .dgp_return(model, X, Xb, n, n.test, sig)
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
    e$method_coverage <- method_coverage
    e$method_post_pred <- method_post_pred
    
    # init cluster
    cl <- parallel::makeCluster(getOption("cl.cores", CORES), outfile="")
    parallel::clusterExport(cl, ls(e, all.names=TRUE), envir=e)
    on.exit(parallel::stopCluster(cl))
    
    # run the method
    res <- parallel::parSapply(cl, 1:runs, function(run) {
	cat(run, "\t")
	d <- dgp(n, p, g, s, bmax, pars, seed=run)
	method(d, method_par)
    })

    return(t(res))
}


m_gsvb <- function(d, m_par=list(family="gaussian", lambda=0.5, a0=1, b0=100, 
	a_t=1e-3, b_t=1e-3, diag_covariance=TRUE, intercept=TRUE))
{
    tryCatch({
	fit.time <- system.time({
	    fit <- gsvb::gsvb.fit(d$y, d$X, d$groups, 
		family=m_par$family, intercept=m_par$intercept, 
		diag_covariance=m_par$diag_covariance, lambda=m_par$lambda,
		a0=m_par$a0, b0=m_par$b0, tau_a0=m_par$a_t, tau_b0=m_par$b_t,
		niter=500, track_elbo=FALSE)
	})

	active_groups <- rep(0, length(unique(d$groups)))
	active_groups[d$active_groups] <- 1
	if (m_par$intercept) {
	    res <- method_summary(d$b, active_groups, fit$beta_hat[-1], fit$g[-1], 0.5)
	} else {
	    res <- method_summary(d$b, active_groups, fit$beta_hat, fit$g, 0.5)
	}
	
	coverage.beta <- method_coverage(d, fit, "gsvb", prob = 0.95)

	if (!is.null(d$test)) {
	    coverage.pp <- method_post_pred(d, fit, method="gsvb", 
		quantiles=c(0.025, 0.975), return_samples=FALSE)

	    return(c(unlist(res), unlist(fit.time[3]), unlist(coverage.beta),
		     unlist(coverage.pp)))
	}

	return(c(unlist(res), unlist(fit.time[3]), unlist(coverage.beta)))

    }, error=function(e) 
    {
	cat("error in run: ", d$seed)
	print(e)

	if (!is.null(d$test)) return(rep(NA, 218))
	return(rep(NA, 217))
    })
}


m_spsl <- function(d, m_par=list(family="gaussian", lambda=0.5, a0=1, b0=100, 
    a_t=1e-3, b_t=1e-3, mcmc_samples=10e3, intercept=TRUE))
{
    fit.time <- system.time({
	fit <- spsl::spsl.fit(d$y, d$X, d$groups, family=m_par$family,
	    intercept=m_par$intercept, lambda=m_par$lambda, a_0=m_par$a0, 
	    b_0=m_par$b0, a_t=m_par$a_t, b_t=m_par$b_t, 
	    mcmc_sample=m_par$mcmc_samples)
    })

    active_groups <- rep(0, length(unique(d$groups)))
    active_groups[d$active_groups] <- 1

    if (m_par$intercept) {
	res <- method_summary(d$b, active_groups, fit$beta_hat[-1], fit$g[-1], 0.5)
    } else {
	res <- method_summary(d$b, active_groups, fit$beta_hat, fit$g, 0.5)
    }

    coverage.beta <- method_coverage(d, fit, "spsl", prob = 0.95)

    if (!is.null(d$test)) {
	coverage.pp <- method_post_pred(d, fit, method="spsl", 
	    quantiles=c(0.025, 0.975), return_samples=FALSE)

	return(c(unlist(res), unlist(fit.time[3]), unlist(coverage.beta),
		 unlist(coverage.pp)))
    }

    return(c(unlist(res), unlist(fit.time[3]), unlist(coverage.beta)))
}


m_ssgl <- function(d, m_par=list(family="gaussian", l0=20, l1=1, a0=1, b0=100)) 
{
    tryCatch({
	fit.time <- system.time({
	    fit <- sparseGAM::SSGL(d$y, d$X, d$X, d$groups, family=m_par$family,
		lambda0=m_par$l0, lambda1=m_par$l1, a=m_par$a0, b=m_par$b0)
	})

	active_groups <- rep(0, length(unique(d$groups)))
	active_groups[d$active_groups] <- 1
	res <- method_summary(d$b, active_groups, fit$beta, fit$classifications, 0.5)
	return(c(unlist(res), unlist(fit.time[3])))
    }, error=function(e) 
    {
	cat("error in run: ", d$seed)
	print(e)

	return(rep(NA, 213))
    })
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
	return_samples=FALSE, samples=1e4) 
{
    if (method == "spsl") {
	pp <- spsl::spsl.predict(fit, d$test$X, quantiles=quantiles,
	    return_samples=return_samples)
    } 
    if (method == "gsvb") {
	pp <- gsvb::gsvb.predict(fit, d$test$X, samples=samples,
	    quantiles=quantiles, return_samples=return_samples)
    }
    coverage <- mean((pp$quantiles[1, ] <= d$test$y) & (d$test$y <= pp$quantiles[2, ]))

    return(list(coverage=coverage))
}


method_coverage <- function(d, fit, method, prob=0.95)
{
    if (method == "spsl") {
	ci <- spsl::spsl.credible_intervals(fit, prob=prob)
    }
    if (method == "gsvb") {
	ci <- gsvb::gsvb.credible_intervals(fit, prob=prob)
    }

    l <- ci[ , 1]
    u <- ci[ , 2]
    dirac <- ci[ , 3]
    
    b <- if (fit$parameters$intercept) c(0, d$b) else d$b
    coverage <- (((l <= b) & (b <= u)) | ((b == 0) & (dirac == T)))
    bs <- b != 0

    return(list(
	coverage.non_zero = mean(coverage[bs]),
	length.non_zero = mean(u[bs] - l[bs]),
	coverage.zero = mean(coverage[!bs]),
	length.zero =  mean(u[!bs] - l[!bs])
    ))
}


# ----------------------------------------
# Misc
# ----------------------------------------
# import an environment variable ex: cores <- d("CORES", 8)
read.env <- function(evar, default) 
{
    if(Sys.getenv(evar) != "") as(Sys.getenv(evar), class(default)) else default
}

