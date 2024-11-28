# -------------------------------------------------------------------------------
# 			Analysis of Rat Eye Data
#
# 	ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5680/matrix/
#
# -------------------------------------------------------------------------------
rm(list=ls())
set.seed(1)

library(AnnotationDbi) 		# BiocManager::install("AnnotationDbi")
library(rat2302.db) 		# BiocManager::install("rat2302.db")
library(splines)
library(data.table)
library(SSGL)
library(foreach)
library(data.table)
# library(doParallel)

source("./00-functions.R")


# -------------------------------------------------------------------------------
# 				Data prep
# -------------------------------------------------------------------------------
dat = data.table::fread("../../data/rat/GSE5680_series_matrix.txt", 
			header=F, fill=T)

col_names = dat[1, -1]
row_names = as.character(dat[-1, 1])
row_names = unlist(dat[-1, 1])

X = dat[-1, -1]
rownames(X) = NULL
colnames(X) = NULL
X = apply(X, 2, as.numeric)

gene_name_table = AnnotationDbi::select(rat2302.db, row_names, 
	c("SYMBOL","ENTREZID", "GENENAME"))
gene_name_table = gene_name_table[!duplicated(gene_name_table$PROBEID), ]

# set up outcome and design matrix
trim32_index = which(row_names == gene_name_table[which(gene_name_table$SYMBOL == "Trim32"), ]$PROBEID)
# trim32_index = which(gene_name_table$SYMBOL == "Trim32")
Y = X[trim32_index, ]
X = X[-trim32_index, ]

# For a probe to be considered expressed, the maximum expression value observed 
# for that probe among the 120 F2 rats was required to be greater than the 
# 25th percentile of the entire set of expression values. 
q25 = quantile(X, 0.25, na.rm=T)
keep = which(apply(X, 1, function(x) max(x) >= q25))
X = X[keep, ]
row_names = row_names[keep]

# keep 1000 genes with highest variance (on the log scale)
X_var = apply(log(X), 1, var)
keep = order(X_var, decreasing = T)[1:5000]
X = X[keep, ]
row_names = row_names[keep]

# get gene symbols from probe ids
gene_name_table = AnnotationDbi::select(rat2302.db, row_names, 
	c("SYMBOL","ENTREZID", "GENENAME"))
gene_name_table = gene_name_table[!duplicated(gene_name_table$PROBEID), ]

X = t(X)

# -------------------------------------------------------------------------------
# 				Splines setup
# -------------------------------------------------------------------------------
create_train_test = function(y, X, tr, ts)
{
    p = ncol(X)
    n_tr = length(tr)
    n_ts = length(ts)

    m = 3 	# basis expansion size
    Z_tr = matrix(NA, nrow=n_tr, ncol=m * p)
    Z_ts = matrix(NA, nrow=n_ts, ncol=m * p)

    for (j in 1:p) {
        basis = splines::ns(X[tr, j], df=m)
        Z_tr[ , m * (j - 1) + 1:m]  = basis
        Z_ts[ , m * (j - 1) + 1:m]  = predict(basis, X[ts, j]) 
    }

    # Package dataset
    d = list(y_train=y[tr], X_train=Z_tr, y_test=y[ts], X_test=Z_ts, 
             groups=rep(1:p, each=m), n=n_tr, M=p, m=m, p=ncol(Z_tr))
}

rm(list=c("q25", "dat", "X_var", "row_names", "col_names",
	  "trim32_index", "keep"))

# -------------------------------------------------------------------------------
# 		   Functions to produce a orthonormal design
#
# -------------------------------------------------------------------------------
std = function(d)
{
    n = nrow(d$X)
    M = length(unique(d$groups))

    center  = colMeans(d$X)
    X. = sweep(d$X, 2, center)
        
    scale = list()

    for (g in 1:M)
    {
	G = which(d$groups == g)

	if (length(G) == 1)
	{
	    scl = sqrt(sum((X.[ , 1])^2) / m)
	    X.[ , G] = X.[ , G] / scl

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
		fit$s[G] = fit$s[G] * d$scaled[[g]]
	    } else {
		fit$s[[g]] = fit$s[[g]] * d$scaled[[g]]^2
	    }
	} else 
	{
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

mse = function(x, y) mean((x - y)^2)


# -------------------------------------------------------------------------------
# 				X-validate
# -------------------------------------------------------------------------------
folds = 10
results = array(NA, dim=c(4, folds, 2))
models = list(diag_cov=list(), full_cov=list())

# train the models 
# you can skip this and load the pre-trianed models by calling:
# 
# load("../../rdata/application/rat/models_15k.RData")
#

for (fold in 1:folds) 
{
    xval = cv(nrow(X), fold, folds=10, random_order=FALSE)
    
    tr = xval$train
    ts = xval$test

    d = create_train_test(Y, X, tr, ts)
    d.train = list(y=d$y_train, X=d$X_train, groups=d$groups)
    d.test  = list(y=d$y_test,  X=d$X_test,  groups=d$groups)
    
    # fit the model
    d.train = std(d.train)

    f1 = gsvb::gsvb.fit(d.train$y., d.train$X., d.train$groups, intercept=FALSE, 
			diag_covariance=TRUE, track_elbo=FALSE)
    f2 = gsvb::gsvb.fit(d.train$y., d.train$X., d.train$groups, intercept=FALSE, 
			diag_covariance=FALSE, track_elbo=FALSE)

    # save the models
    models$diag_cov[[fold]] = f1
    models$full_cov[[fold]] = f2 

    # get model
    f1 = rescale_fit(models$diag_cov[[fold]], d.train)
    f2 = rescale_fit(models$full_cov[[fold]], d.train)
    
    f1$intercept = mean(d.test$y - d.test$X %*% (f1$mu * (f1$g[d.train$groups] > 0.5)))
    f2$intercept = mean(d.test$y - d.test$X %*% (f2$mu * (f2$g[d.train$groups] > 0.5)))

    # evaluate the model
    results[1, fold, 1] = mse(d.test$y, 
	    d.test$X %*% (f1$mu * (f1$g[d.train$groups] > 0.5)) + f1$intercept)

    results[1, fold, 2] = mse(d.test$y, 
	    d.test$X %*% (f2$mu * (f2$g[d.train$groups] > 0.5)) + f2$intercept)

    results[2, fold, 1] = sum(f1$g > 0.5)
    results[2, fold, 2] = sum(f2$g > 0.5)

    # compute the PP coverage
    f1.pred = gsvb::gsvb.predict(f1, d.test$X)

    results[3, fold, 1] = 
        mean(d.test$y - f1$intercept >= f1.pred$quantiles[1, ] & 
         d.test$y - f1$intercept <= f1.pred$quantiles[2, ])

    results[4, fold, 1] = 
        mean(abs(f1.pred$quantiles[2, ] - f1.pred$quantiles[1, ] ))


    f2.pred = gsvb::gsvb.predict(f2, d.test$X)

    results[3, fold, 2] = 
        mean(d.test$y - f2$intercept >= f2.pred$quantiles[1, ] & 
         d.test$y - f2$intercept <= f2.pred$quantiles[2, ])

    results[4, fold, 2] = 
        mean(abs(f2.pred$quantiles[2, ] - f2.pred$quantiles[1, ] ))
}


# Note: we use the exact same mechanism of perfroming X-validation as SSGL
# the following takes a while to run ~ 12 hours
# Cross-validation loop to find the optimal value of lambda0 for SSGL
lambda0_values <- seq(0.1, 10, by=0.5)
cv_results <- matrix(NA, nrow=length(lambda0_values), ncol=10)

for (i in seq_along(lambda0_values)) {
    lambda0 <- lambda0_values[i]
    
    for (fold in 1:10) {
        xval <- cv(nrow(X), fold, folds=10, random_order=FALSE)
        tr <- xval$train
        ts <- xval$test

        d <- create_train_test(Y, X, tr, ts)
        d.train <- list(y=d$y_train, X=d$X_train, groups=d$groups)
        d.test  <- list(y=d$y_test, X=d$X_test, groups=d$groups)

        # standardization happens inside SSGL
        f.s <- SSGL::SSGL(d.train$y, d.train$X, d.train$groups, family="gaussian",
                          d.test$X, lambda0=lambda0)

        cv_results[i, fold] <- mse(d.test$y, f.s$Y_pred)
    }
}

# Calculate the mean MSE for each lambda0 value
mean_mse <- rowMeans(cv_results)
which.min(mean_mse)

# Find the optimal lambda0 value
optimal_lambda0 <- lambda0_values[which.min(mean_mse)]
cat("Optimal lambda0 value: ", optimal_lambda0, "\n")

# this was the optimal value of lambda0 from xval
lambda0 = 1.6

# extract the results for the mse minimizing lambda0
# note: ssglcv does not provide the models, just the param
# so to get the model size and which coefficients are important
# we need to run it again
ssgl_models = list()
ssgl_results = matrix(NA, nrow=2, ncol=10)

for (fold in 1:10) 
{
    xval = cv(nrow(X), fold, folds=10, random_order=FALSE)
    tr = xval$train
    ts = xval$test

    d = create_train_test(Y, X, tr, ts)
    d.train = list(y=d$y_train, X=d$X_train, groups=d$groups)
    d.test  = list(y=d$y_test, X=d$X_test, groups=d$groups)

    # standardization happens inside SSGL
    f.s = SSGL::SSGL(d.train$y, d.train$X, d.train$groups, family="gaussian",
		     d.test$X, lambda0=lambda0)

    ssgl_models[[fold]] = f.s

    f.s = ssgl_models[[fold]]

    ssgl_results[1, fold] = mse(d.test$y, f.s$Y_pred)
    ssgl_results[2, fold] = sum(f.s$classifications)
}

# rowMeans(ssgl_results)

# save the models
save(list=c("models", "ssgl_models"), 
     file="../../rdata/application/rat/models_15k.RData")

# save the results
save(list=c("results", "ssgl_results"), 
     file="../../rdata/application/rat/results.RData")



# -------------------------------------------------------------------------------
# 				Print the results
# -------------------------------------------------------------------------------
load("../../rdata/application/rat/results.RData")
load("../../rdata/application/rat/models_15k.RData")

cat(apply(results, c(1, 3), function(x) 
	sprintf("%.4f (%.3f)", mean(x[ ! (is.na(x) | is.nan(x)) ]), 
			       sd(x[ ! (is.na(x) | is.nan(x)) ]))
))
cat(apply(ssgl_results, 1, function(x) 
	sprintf("%.4f (%.3f)", mean(x[ ! (is.na(x) | is.nan(x)) ]), 
			       sd(x[ ! (is.na(x) | is.nan(x)) ]))
))

cbind(gene_name_table[
as.numeric(names(table(unlist(lapply(models$diag_cov, function(f) which(f$g > 0.5))))))
, c(1,2)],
table(unlist(lapply(models$diag_cov, function(f) which(f$g > 0.5)))))

cbind(gene_name_table[
as.numeric(names(table(unlist(lapply(models$full_cov, function(f) which(f$g > 0.5))))))
, c(1, 2)],
table(unlist(lapply(models$full_cov, function(f) which(f$g > 0.5)))))

cbind(gene_name_table[
as.numeric(names(table(unlist(lapply(ssgl_models,function(f) which(f$classifications != 0)))))), c(1,2) ],
table(unlist(lapply(ssgl_models, function(f) which(f$classifications != 0)))))



# -------------------------------------------------------------------------------
# 			Fit the models to the full dataset
# -------------------------------------------------------------------------------
# standardize the data
d = std(d)

f.ssgl = SSGL(d$y, d$X, d$X, 1, 100, d$groups)
f1 = gsvb::gsvb.fit(d$y., d$X., d$groups, intercept=FALSE, diag_covariance=TRUE)
f2 = gsvb::gsvb.fit(d$y., d$X., d$groups, intercept=FALSE, diag_covariance=FALSE)

save(list=c("f1", "f2"), file="../../rdata/application/rat/models_15k_full_gsvb.RData")

gene_name_table[f1$g > 0.5, ]
gene_name_table[f2$g > 0.5, ]

mse(d$y., d$X. %*% f1$beta_hat)
mse(d$y., d$X. %*% f2$beta_hat)

f1.pred = gsvb::gsvb.predict(f1, d$X.)
mean(d$y. >= f1.pred$quantiles[1, ] & d$y. <= f1.pred$quantiles[2, ])
mean(abs(f1.pred$quantiles[2, ] - f1.pred$quantiles[1, ] ))

f2.pred = gsvb::gsvb.predict(f2, d$X.)
mean(d$y. >= f2.pred$quantiles[1, ] & d$y. <= f2.pred$quantiles[2, ])
mean(abs(f2.pred$quantiles[2, ] - f2.pred$quantiles[1, ] ))
    
