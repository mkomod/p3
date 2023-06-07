# -------------------------------------------------------------------------------
# 			Analysis of Rat Eye Data
#
# 	ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5680/matrix/
#
# -------------------------------------------------------------------------------

library(AnnotationDbi) 		# BiocManager::install("AnnotationDbi")
library(rat2302.db) 		# BiocManager::install("rat2302.db")
library(splines)
library(data.table)
library(SSGL)


rm(list=ls())
set.seed(1)

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

# get gene symbols from probe ids
gene_name_table = AnnotationDbi::select(rat2302.db, row_names, 
	c("SYMBOL","ENTREZID", "GENENAME"))
gene_name_table = gene_name_table[!duplicated(gene_name_table$PROBEID), ]

# set up outcome and design matrix
trim32_index = which(gene_name_table$SYMBOL == "Trim32")
Y = X[trim32_index, ]
X = X[-trim32_index, ]

# For a probe to be considered expressed, the maximum expression value observed 
# for that probe among the 120 F2 rats was required to be greater than the 
# 25th percentile of the entire set of expression values. 
q25 = quantile(X, 0.25, na.rm=T)
keep = which(apply(X, 1, function(x) max(x) >= q25))
X = X[keep, ]
gene_name_table = gene_name_table[keep, ]

# keep 1000 genes with highest variance (on the log scale)
X_var = apply(log(X), 1, var)
keep = order(X_var, decreasing = T)[1:1000]
X = X[keep, ]
gene_name_table = gene_name_table[keep, ]



# -------------------------------------------------------------------------------
# 				Splines setup
# -------------------------------------------------------------------------------
X = t(X)
X = scale(X)
p = ncol(X)
n = nrow(X)

m = 3 	# basis expansion size
Z = matrix(nrow=n, ncol=0)

for (j in 1:p)
    Z = cbind(Z, ns(X[, j], df=m))

Z = scale(Z, center=T, scale=F)



# -------------------------------------------------------------------------------
# 				Fit the models
# -------------------------------------------------------------------------------
groups = rep(1:p, each=m)
y = scale(Y, center=TRUE, scale=FALSE)

# Fit GSVB
f1 = gsvb::gsvb.fit(y, Z, groups, intercept=FALSE)
f2 = gsvb::gsvb.fit(y, Z, groups, intercept=FALSE, diag_covariance=FALSE)

plot(f1$g)
plot(f1$beta)
gene_name_table[ which(f1$g > 0.5), ]

plot(f2$g)
plot(f2$beta)
gene_name_table[ which(f2$g > 0.5), ]


# Fit SSGL
lambda1 = 1
lambda0seq = seq(1, 500, length.out = 100)

modSSGLcv = SSGLcv(Y=y, X=Z, lambda1=lambda1, lambda0seq=lambda0seq,
    groups=groups, a=1, b=p, nFolds=10, M=10, error=0.001)
which.min(modSSGLcv$CVerror)



# -------------------------------------------------------------------------------
# 				X-validate
# -------------------------------------------------------------------------------
# mod.gsvb = list()
# mod.ssgl = list()
# load("../../rdata/application/rat/cv_results.RData")

cv.error = matrix(0, nrow=2, ncol=10)
mod.size = matrix(0, nrow=2, ncol=10) 
pp.coverage = matrix(0, nrow=1, ncol=10)

for (fold in 1:10) 
{
    tt = cv(nrow(X), fold, folds=10)
    tr = tt$train

    # f.gsvb = gsvb::gsvb.fit(y[tr], Z[tr, ], groups,
			# intercept=FALSE, diag_covariance=FALSE, 
			# track_elbo=FALSE)
    # mod.gsvb[[fold]] = f.gsvb

    # f.ssgl = SSGL(Y=y[tr], X=Z[tr, ], lambda1=1, lambda0=91,
		   # groups = groups,a=1, b=p, M=10, error=0.001,
		   # theta = 0.5)
    # mod.ssgl[[fold]] = f.ssgl

    f.gsvb = mod.gsvb[[fold]]
    f.ssgl = mod.ssgl[[fold]]

    
    ts = tt$test
    cv.error[1, fold] = mean((y[ts] - Z[ts, ] %*% f.gsvb$beta_hat)^2)
    cv.error[2, fold] = mean((y[ts] - f.ssgl$intercept - Z[ts, ] %*% f.ssgl$beta)^2)

    mod.size[1, fold] = sum(f.gsvb$g > 0.5)
    mod.size[2, fold] = sum(f.ssgl$beta != 0) / m
    
    gsvb.pred = gsvb::gsvb.predict(f.gsvb, Z[ts, ])
    pp.coverage[1, fold] = mean(y[ts] >= gsvb.pred$quantiles[1, ] & 
				y[ts] <= gsvb.pred$quantiles[2, ])
}

# save(list=c("mod.gsvb", "mod.ssgl", "cv.error", "mod.size"),
#      file="../../rdata/application/rat/cv_results.RData")

rowMeans(cv.error)
apply(cv.error, 1, sd)

rowMeans(mod.size)
apply(mod.size, 1, sd)

rowMeans(pp.coverage)
sd(pp.coverage)


# -------------------------------------------------------------------------------
# 				Tables
# -------------------------------------------------------------------------------
for (j in 1:2)
    cat(sprintf("%.4f (%.3f) & %.2f (%.3f) \\\\ \n",
    mean(cv.error[j, ]), sd(cv.error[j, ]), 
    mean(mod.size[j, ]),sd(mod.size[j, ])))



# -------------------------------------------------------------------------------
# 				Figures
# -------------------------------------------------------------------------------
BETA = gsvb::gsvb.sample(f2)
layout(matrix(1:4, nrow=2))

for (j in which(f2$g > 0.5)) {
    o = order(X[ , j])
    
    mc_sample = Z[o , groups == j] %*% BETA$beta[groups == j, ]
    quant = apply(mc_sample, 1, quantile, probs=c(0.005, 0.995))

    plot(X[o , j], Z[o , groups == j] %*% f2$beta_hat[groups == j], type="b", lwd=3, 
	 ylim=range(quant), ylab=latex2exp::TeX(paste0("$f_", j, "(x)$")), xlab="x", 
	 pch=4)
    lines(X[o, j], quant[1, ], lty=2)
    lines(X[o, j], quant[2, ], lty=2)
    grid(10, 10)
}
plot.new()


pred = gsvb::gsvb.predict(f2, Z)


plot(y, pred$mean)
y >= pred$quantiles[1, ] & y <= pred$quantiles[2, ]
mean(y >= pred$quantiles[1, ] & y <= pred$quantiles[2, ])

