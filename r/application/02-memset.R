# -------------------------------------------------------------------------------
# 			Analysis of MEMSET data
#
#		http://hollywood.mit.edu/burgelab/maxent/ssdata/
#
# -------------------------------------------------------------------------------
rm(list=ls())
set.seed(1)

library(ROSE)
library(gglasso)
library(sparseGAM)
library(gsvb)

source("./00-functions.R")


# -------------------------------------------------------------------------------
# 			Load and process data
# -------------------------------------------------------------------------------
dat.train = read.delim("../../data/memset/trainfile5_hs", header=F)

n.train = nrow(dat.train) / 2
train.X = matrix(nrow=n.train, ncol=7)
train.Y = matrix(nrow=n.train, ncol=1)

for (i in 1:n.train) {
    train.X[i, ] = strsplit(dat.train[i*2, ], "")[[1]]
    train.Y[i, ] = stringi::stri_sub(dat.train[(i-1)*2 + 1, ], 3)
}

# code (0, 5) -> (0, 1)
train.Y = matrix(as.numeric(train.Y))
train.Y[train.Y == 0, ] = 0
train.Y[train.Y == 5, ] = 1

# package up data into dataframe
d.train = data.frame(y=train.Y, train.X)

# Sub-sample the data
set.seed(1234)
keep0 = sort(sample(which(d.train$y == 0), 1500, replace=FALSE))
keep1 = sort(sample(which(d.train$y == 1), 1500, replace=FALSE))

# create a balanced dataset
d.balanced = d.train[c(keep0, keep1), ]

# create the model matrix
X = model.matrix(y ~ .^3, d.balanced)
y = d.balanced$y
yy = y
yy[yy == 0] = -1

# create the groups
cnames = colnames(X)
cnames = stringi::stri_replace_all(cnames, "", regex="[a|g|c|t|X]")
cnames[1] = "Intercept"
ucnames = unique(cnames)

groups = double(length(cnames))
for (i in 1:length(ucnames))
    groups[which(cnames == ucnames[i])] = i


# package up the data
# yy: is a version of the response with {-1, 1} labels
# XX: is a rescaled version of the design matrix that  plays nice with 
#     gsvb when the covariance is block diag
d = list(y=y, yy=yy, X=X, XX=X*5, groups=groups)


rm(list=c("d.train", "keep1", "keep0", "i", "y", "X", "yy", 
	  "train.X", "train.Y", "n.train", "cnames"))


d = std(d, FALSE)
d$O = grpreg:::orthogonalize(d$X, d$groups)


# -------------------------------------------------------------------------------
# 				Fit models
# -------------------------------------------------------------------------------
# run this to skip all the model fitting
load("../../rdata/application/memset/models.RData")


# -----------------------------------------------------------------
#   Fit group LASSO
cv.glasso = gglasso::cv.gglasso(d$X, d$yy, d$groups, loss="logit",
	    pred.loss="misclass", intercept=FALSE)
lambda.min = cv.glasso$lambda[which.min(cv.glasso$cvlo)] # lambda.min =  0.001766377
plot(cv.glasso)

f.gg = gglasso::gglasso(d$X, d$yy, loss="logit", lambda=lambda.min, intercept=FALSE)

plot(f.gg$beta)
plot.norms(f.gg$beta, d$groups, ucnames)
sum(nzero.groups(f.gg$beta, groups))


# -----------------------------------------------------------------
#   Fit group spike-and-slab LASSO
cv.ssgl = sparseGAM::cv.SSGL(d$y, d$X, d$groups, family="binomial", 
	    nlambda0=10)
f.ssgl.full = sparseGAM::SSGL(d$y, d$X, d$X, groups, family="binomial")

matplot(t(f.ssgl.full$beta), type="l")

# fit ssgl
lambda0.min = 4
f.ssgl = sparseGAM::SSGL(d$y, d$X, d$X[1:10, ], groups, family="binomial",
		lambda0=lambda0.min)
f.ssgl$beta[1] = f.ssgl$beta0
f.ssgl$class[1] = 1

# -----------------------------------------------------------------
#   Fit GSVB
f1 = gsvb::gsvb.fit(d$y, d$X., d$groups, family="binomial-jaakkola", lambda=1,
	    diag_covariance=TRUE, intercept=FALSE, track_elbo=FALSE, niter=1000)

f1. = rescale_fit(f1, d)

f2 = gsvb::gsvb.fit(d$y, d$X., d$groups, family="binomial-jaakkola", lambda=1,
	    diag_covariance=FALSE, intercept=FALSE, track_elbo=FALSE, niter=1500)

f2. = rescale_fit(f2, d)

# -----------------------------------------------------------------
#   Save
save(list=c("cv.glasso", "f.gg", 
	    "cv.ssgl", "f.ssgl.full", 
	    "f.ssgl",
	    "f1", "f2"),
     file="../../rdata/application/memset/models.RData")

# -----------------------------------------------------------------------------
# 			Plot the norms of the methods
# -----------------------------------------------------------------------------

q1. = apply(comp.norms(f1., d$groups), 1, quantile, probs=c(0.025, 0.975))
q2. = apply(comp.norms(f2., d$groups), 1, quantile, probs=c(0.025, 0.975))
n.g = nzero.groups(f.gg$beta, d$groups)



pdf("../figs/memset_norms.pdf", 9, 6)
layout(matrix(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,4,4), ncol=1))

# GSVB-D
par(mar=c(0, 4, 0.1, 1), family="Times")
plot.norms(f1.$beta, groups, ylim=c(0, 18.0), pch=20, xaxt="n", las=1,
	   col=ifelse(f1.$g > 0.5, 1, "lightgrey"), ylab="GSVB-D")
arrows((1:max(groups)) * (-1)^((f1.$g > 0.5) - 1), q1.[1, ], y1=q1.[2, ], 
       code=0, col=1, lwd=2)
selected = (f1.$g > 0.5 & f2.$g > 0.5) | (f1.$g > 0.5 & f.ssgl$class) | (f1.$g > 0.5 & n.g)
abline(v = which(selected), lty=2, lwd=0.5, col=2)
grid(nx=0, ny=5)

# GSVB-B
par(mar=c(0, 4, 0, 1))
plot.norms(f2.$beta, groups, ylim=c(0, 230), pch=20, xaxt="n", las=1,
	   col=ifelse(f2.$g > 0.5, 1, "lightgrey"), ylab="GSVB-B")
arrows((1:max(groups)) * (-1)^((f2.$g > 0.5) - 1), 
       q2.[1, ], y1=q2.[2, ], code=0, col=1, lwd=2)
selected = (f2.$g > 0.5 & f1.$g > 0.5) | (f2.$g > 0.5 & f.ssgl$class) | (f2.$g > 0.5 & n.g)
abline(v = which(selected), lty=2, lwd=0.5, col=2)
grid(nx=0, ny=5)

# SSGL
par(mar=c(0, 4, 0, 1))
plot.norms(f.ssgl$beta, groups, ylim=c(0, 18), pch=20, xaxt="n", las=1,
	   col=ifelse(f.ssgl$class, 1, "lightgrey"), ylab="SSGL")
selected = (f2.$g > 0.5 & f.ssgl$class) | (f1.$g > 0.5 & f.ssgl$class) | (n.g & f.ssgl$class)
abline(v = which(selected), lty=2, lwd=0.5, col=2)
grid(nx=0, ny=5)

# Group LASSO
par(mar=c(3, 4, 0, 1))
plot.norms(f.gg$beta, groups, ylim=c(0, 13), pch=20, las=1, xaxt="n",
	   col=ifelse(n.g, 1, "lightgrey"), ylab="Group LASSO")
selected = (n.g & f1.$g > 0.5) | (f2.$g > 0.5 & n.g) | (f.ssgl$class & n.g)
abline(v = which(selected), lty=2, lwd=0.5, col=2)
grid(nx=0, ny=5)

axis(1, 1:max(groups), F, lwd=0, lwd.ticks=0.8)
text(x=1:max(groups), par("usr")[3] - 0.8,
     labels = ucnames, srt = -45, 
     pos = 1, xpd = TRUE, cex=0.82)


dev.off()



# -------------------------------------------------------------------------------
# 			Create the test data
# -------------------------------------------------------------------------------
dat.test = read.delim("../../data/memset/testfile5_hs", header=F)

n.test = nrow(dat.test) / 2
test.X = matrix(nrow=n.test, ncol=7)
test.Y = matrix(nrow=n.test, ncol=1)

for (i in 1:n.test) {
    test.X[i, ] = strsplit(dat.test[i*2, ], "")[[1]]
    test.Y[i, ] = stringi::stri_sub(dat.test[(i-1)*2 + 1, ], 3)
}

test.Y = matrix(as.numeric(test.Y))
test.Y[test.Y == 0, ] = 0
test.Y[test.Y == 5, ] = 1

df.test = data.frame(y=test.Y, test.X)

# Sub-sample the data
set.seed(1234)
keep0 = sort(sample(which(df.test$y == 0), 4208, replace=FALSE))
keep1 = sort(sample(which(df.test$y == 1), 4208, replace=FALSE))
# keep0 = sort(sample(which(df.test$y == 0), 1500, replace=FALSE))
# keep1 = sort(sample(which(df.test$y == 1), 1500, replace=FALSE))

d.balanced = df.test[c(keep0, keep1), ]

# create the model matrix
X = model.matrix(y ~ .^3, d.balanced)
y = d.balanced$y
yy = y
yy[yy == 0] = -1


# create the groups
cnames = colnames(X)
cnames = stringi::stri_replace_all(cnames, "", regex="[a|g|c|t|X]")
cnames[1] = "Intercept"
ucnames = unique(cnames)

groups = double(length(cnames))
for (i in 1:length(ucnames))
    groups[which(cnames == ucnames[i])] = i


# package up the data
d.test = list(y=y, yy=yy, X=X, XX=X*5, groups=groups)


# -------------------------------------------------------------------------------
# 			Compare the models on test data
# -------------------------------------------------------------------------------
dd.test = d.test
results = array(NA, dim=c(6, 10, 4))

for (fold in 1:10) 
{
    xval = cv(4208, fold, 10, random_order=FALSE)
    idx = c(xval$test, 4208 + xval$test)
    d.test = list(y = dd.test$y[ idx] , X = dd.test$X[idx, ])

    glas.pred = 1 / (1 + exp(- d.test$X %*% f.gg$beta))
    ssgl.pred = 1 / (1 + exp(- d.test$X %*% f.ssgl$beta))
    gsvb.d.pred = 1 / (1 + exp(- d.test$X %*% (f1.$mu * (f1.$g[groups] > 0.5))))
    gsvb.b.pred = 1 / (1 + exp(- d.test$X %*% (f2.$mu * (f2.$g[groups] > 0.5))))

    results[ , fold, 1] = cat.res(d.test$y, as.numeric(gsvb.d.pred), 
			       f1.$mu * (f1.$g[groups] > 0.5), groups)
    results[ , fold, 2] = cat.res(d.test$y, as.numeric(gsvb.b.pred), 
			       f2.$mu * (f2.$g[groups] > 0.5), groups)
    results[ , fold, 3] = cat.res(d.test$y, as.numeric(ssgl.pred), f.ssgl$beta, groups)
    results[ , fold, 4] = cat.res(d.test$y, as.numeric(glas.pred), f.gg$beta, groups)

}

apply(results, c(1, 3), function(x) {
	for (i in x)
	    sprintf("%.3f (%.3f)", mean(x), sd(x)) 
})
t(round(apply(results, c(1,3), mean), 3))
t(round(apply(results, c(1,3), sd), 4))

max(sapply(taus, function(tau) cor(d.test$y, glas.pred > tau)))
max(sapply(taus, function(tau) cor(d.test$y, ssgl.pred > tau)))
max(sapply(taus, function(tau) cor(d.test$y, gsvb.d.pred > tau))) 
max(sapply(taus, function(tau) cor(d.test$y, gsvb.b.pred > tau)))


# ROSE::roc.curve(d.test$y, as.numeric(glas.pred), col=3, add=T)
ROSE::roc.curve(d.test$y, as.numeric(glas.pred))
ROSE::roc.curve(d.test$y, as.numeric(ssgl.pred))
ROSE::roc.curve(d.test$y, as.numeric(gsvb.d.pred), col=2, add=T)
ROSE::roc.curve(d.test$y, as.numeric(gsvb.b.pred), col=3, add=T)

{
cat("GSVB-B &"); 
cat.res(d.test$y, as.numeric(gsvb.d.pred), f1.$mu * (f1.$g[groups] > 0.5), groups)
cat("GSVB-D &"); 
cat.res(d.test$y, as.numeric(gsvb.b.pred), f2.$mu * (f2.$g[groups] > 0.5), groups)
cat("SSGL &"); 
cat.res(d.test$y, as.numeric(ssgl.pred), f.ssgl$beta, groups)
cat("Group LASSO &"); 
cat.res(d.test$y, as.numeric(glas.pred), f.gg$beta, groups)
}

# correlation between true and predicted for different threshold

ucnames[f1.$g > 0.5]
ucnames[f2.$g > 0.5]
