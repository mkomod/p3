

cv = function (n, fold, folds = 10, random_order = TRUE, seed = 1)
{
    # perfrom cross validation

    if (fold < 0 || fold > folds)
        stop("fold must be between 0 and folds")
    if (folds > n)
        stop("Cannot have more folds than n")
    if (random_order) {
        set.seed(seed)
        index <- sample(1:n)
    }
    else {
        index <- 1:n
    }
    k <- floor(n/folds)
    test <- index[((fold - 1) * k + 1):(fold * k)]
    train <- setdiff(index, test)
    return(list(test = test, train = train))
}


comp.norms = function(fit, groups, probs=c(0.025, 0.0975)) 
{
    BETA = gsvb::gsvb.sample(fit, samples = 5000)$beta

    print(dim(BETA) )

    NORMS = apply(BETA, 2, function(beta) {
	norms = c()
	for (g in unique(groups))
	    norms = c(norms, norm(matrix(beta[groups == g])))
	norms
    })

    return(NORMS)
}


plot.norms = function(beta, groups, group.names=NULL, ...) 
{
    norms = c()
    for (g in unique(groups))
	norms = c(norms, norm(matrix(beta[groups == g])))
    
    if (!is.null(group.names)) {
	plot(norms, xaxt="n", ...)
	axis(1, at=1:64, labels=group.names)
    } else {
	plot(norms, ...)
    }
    norms
}



nzero.groups = function(beta, groups) 
{
    sapply(unique(groups), function(g) any(beta[groups == g] != 0))
}


cat.res = function(true, pred, beta, groups)
{
    res = ROSE::accuracy.meas(true, pred)
    auc = ROSE::roc.curve(true, pred, plotit=F)
    msize = sum(nzero.groups(beta, groups))

    cat(sprintf("%.3f & %.3f & %.3f & %.3f & %d \\\\\n",
	res$precision, res$recall, 2 *res$F, auc$auc, msize))

    invisible()
}


std = function(d, centerX=T)
{
    n = nrow(d$X)
    M = length(unique(d$groups))
    
    center  = colMeans(d$X)
    if (centerX) {
	X. = sweep(d$X, 2, center)
    } else {
	X. = d$X
    }
        
    scale = list()

    for (g in 1:M)
    {
	G = which(d$groups == g)

	if (length(G) == 1)
	{
	    scl = sqrt(sum((X.[ , 1])^2) / n)
	    X.[ , G] = X.[ , G] / scl

	    if (scl == 0) {
		X.[ , G] = 1
		scl = 1
	    }

	    scale[[g]] = scl
	} else 
	{
	    SVD = svd((1/n) * crossprod(X.[ , G]))
	    scl = SVD$u %*% diag(1/sqrt(SVD$d))
	    X.[, G] = X.[, G] %*% scl

	    scale[[g]] = scl
	}
    }

    # scale y
    intercept = mean(d$y)
    y. = d$y - intercept

    return(c(d, list(X.=X., y.=y., intercept=intercept, scale=scale, center=center)))
}


rmstd = function(d)
{
    d$X. = d$scale = d$intercept = d$y. = d$intercept = NULL
    d
}


rescale_fit = function(fit, d) 
{
    if (fit$parameters$intercept) stop("not implemented")

    for (g in unique(fit$parameters$groups)) 
    {
	G = which(g == fit$parameters$groups)

	if (length(G) == 1) 
	{
	    fit$mu[G] = fit$mu[G] * d$scale[[g]]
	    if (fit$parameters$diag_covariance) {
		fit$s[G] = fit$s[G] * d$scale[[g]]
	    } else {
		fit$s[[g]] = fit$s[[g]] * d$scale[[g]]^2
	    }
	} else 
	{
	    # fit$mu[G] = d$scale[[g]] %*% fit$mu[G]
	    fit$mu[G] = d$scale[[g]] %*% fit$mu[G]
	    if (fit$parameters$diag_covariance) {
		fit$s[G] = 
		    sqrt( diag( d$scale[[g]] %*% diag(fit$s[G]^2) %*% t(d$scale[[g]]) ))
	    } else {
		fit$s[[g]] = d$scale[[g]] %*% fit$s[[g]] %*% t(d$scale[[g]])
	    }
	}
    }
    fit$beta_hat = fit$mu * fit$g[fit$parameters$groups]
    fit$intercept = mean(d$y - d$X %*% fit$beta_hat)

    return(fit)
}


rescale_beta = function(fit, d)
{
    b = fit$beta_hat

    # for (g in which(fit$g > 0.5))
    for (g in which(fit$g > 0.5))
    {
	G = which(fit$parameters$groups == g)

	if (length(G) == 1) 
	{
	    b[G] = b[G] * d$scale[[g + fit$parameters$intercept]]
	} else 
	{
	    b[G] = d$scale[[g + fit$parameters$intercept]] %*% b[G]
	}
    }

    return(c(fit, list(beta_hat_scaled = b)))
}


rescale_BETA = function(BETA, d) 
{

    for (g in unique(d$groups))
    {
	G = which(d$groups == g)

	if (length(G) == 1) 
	{
	    BETA[G, ] = BETA[G, ] * d$scale[[g]]
	} else 
	{
	    BETA[G , ] = d$scale[[g]] %*% BETA[G, ]
	}
    }

    return(BETA)
}


rescale_intercept = function(fit, d) 
{
    if (!is.null(fit$beta_hat_scaled)) {
	fit$intercept_scaled = mean(d$y - d$X %*% fit$beta_hat_scaled)
    } else {
	stop("rescale beta first")
    }

    return(fit)
}
