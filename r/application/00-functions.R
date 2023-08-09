

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
