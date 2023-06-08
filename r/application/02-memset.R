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

train.Y = matrix(as.numeric(train.Y))
train.Y[train.Y == 0, ] = 0
train.Y[train.Y == 5, ] = 1

d.train = data.frame(y=train.Y, train.X)


# Sub-sample the data
set.seed(1234)
keep0 = sort(sample(which(d.train$y == 0), 1500, replace=FALSE))
keep1 = sort(sample(which(d.train$y == 1), 1500, replace=FALSE))

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



# -------------------------------------------------------------------------------
# 				Fit models
# -------------------------------------------------------------------------------
#   Fit group LASSO
cv.glasso = gglasso::cv.gglasso(d$X, d$yy, d$groups, loss="logit",
	    pred.loss="misclass", intercept=FALSE)
plot(cv.glasso)
lambda.min = cv.glasso$lambda[which.min(cv.glasso$cvlo)] # 0.002865525

f.gg = gglasso::gglasso(d$X, d$yy, loss="logit", lambda=lambda.min, 
	    intercept=FALSE)
plot(f.gg$beta)
plot.norms(f.gg$beta, d$groups, ucnames)
sum(nzero.groups(f.gg$beta, groups))


# -----------------------------------------------------------------
#   Fit group spike-and-slab LASSO
cv.ssgl = sparseGAM::cv.SSGL(d$y, d$X, d$groups, family="binomial", 
	    nlambda0=10)

plot(cv.ssgl$lambda0, cv.ssgl$cve)
plot(cv.ssgl$lambda0, cv.ssgl$cvse)

f.ssgl.full = sparseGAM::SSGL(d$y, d$X, d$X, groups, family="binomial")

f.ssgl = sparseGAM::SSGL(d$y, d$X, d$X[1:10, ], groups, family="binomial",
		lambda0=cv.ssgl$lambda0.min)

# f.ssgl = sparseGAM::SSGL(y, X, X[1:10, ], groups, family="binomial",
# 		lambda0=25)


# -----------------------------------------------------------------
#   Fit GSVB
f.gsvb.d = gsvb::gsvb.fit(d$y, d$X, d$groups, family="binomial-jaakkola",
	    diag_covariance=TRUE, intercept=FALSE, track_elbo=FALSE, niter=1000)

f.gsvb.b = gsvb::gsvb.fit(d$y, d$XX, d$groups, family="binomial-jaakkola",
	    diag_covariance=FALSE, intercept=FALSE, track_elbo=FALSE, niter=1000)

plot(f.gsvb.d$beta)
plot.plot.norms(f.gsvb.d)
sum(f.gsvb.d$g > 0.5)

plot(f.gsvb.b$beta)
plot.norms(f.gsvb.b$beta, groups)
sum(f.gsvb.d$g > 0.5)


# -----------------------------------------------------------------
#   Save everything
save(list=c("f.gsvb.b", "f.gsvb.d", "f.gg", "cv.glasso", "cv.ssgl"),
     file="../../rdata/application/memset/models.RData")


# -------------------------------------------------------------------------------
# 				Compare the models
# -------------------------------------------------------------------------------
ssgl.pred = 1 / (1 + exp(- cbind(1, X) %*% c(f.ssgl$beta0[1], f.ssgl$beta)))
glas.pred = 1 / (1 + exp(- X %*% f.gg$beta))
gsvb.d.pred = gsvb::gsvb.predict(f.gsvb.d, d.test$X)
gsvb.b.pred = gsvb::gsvb.predict(f.gsvb.b, d.test$XX)

table(d$y, glas.pred > 0.5)
table(d$y, ssgl.pred > 0.5)
table(d$y, gsvb.d.pred$mean > 0.5)
table(d$y, gsvb.b.pred$mean > 0.5)



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
keep0 = sort(sample(which(df.test$y == 0), 1500, replace=FALSE))
keep1 = sort(sample(which(df.test$y == 1), 1500, replace=FALSE))

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
glas.pred = 1 / (1 + exp(- d.test$X %*% f.gg$beta))
ssgl.pred = 1 / (1 + exp(- cbind(1, d.test$X) %*% c(f.ssgl$beta0[1], f.ssgl$beta)))
gsvb.d.pred = gsvb::gsvb.predict(f.gsvb.d, d.test$X)
gsvb.b.pred = gsvb::gsvb.predict(f.gsvb.b, d.test$XX)

# table(d.test$y, glas.pred > 0.5)
table(d.test$y, ssgl.pred > 0.5)
table(d.test$y, gsvb.d.pred$mean > 0.5)
table(d.test$y, gsvb.b.pred$mean > 0.5)

# ROSE::roc.curve(d.test$y, as.numeric(glas.pred), col=3, add=T)
ROSE::roc.curve(d.test$y, as.numeric(ssgl.pred))
ROSE::roc.curve(d.test$y, as.numeric(gsvb.d.pred$mean), col=2, add=T)
ROSE::roc.curve(d.test$y, as.numeric(gsvb.b.pred$mean), col=3, add=T)

ROSE::accuracy.meas(d.test$y, as.numeric(ssgl.pred))
ROSE::accuracy.meas(d.test$y, as.numeric(gsvb.d.pred$mean))
ROSE::accuracy.meas(d.test$y, as.numeric(gsvb.b.pred$mean))

ssgl.beta = f.ssgl$beta
ssgl.beta[1] = f.ssgl$beta0
sum(nzero.groups(ssgl.beta, groups) )
sum(nzero.groups(f.gg$beta, groups) )
sum(f.gsvb.d$g > 0.5)
sum(f.gsvb.b$g > 0.5)

plot(f.gsvb.d$g)
points(f.gsvb.b$g, add=T, pch=20)

# plot the l2 norms for each group
plot.norms(ssgl.beta, groups, ucnames)
plot.norms(f.gsvb.d$beta, groups, ucnames)
plot.norms(f.gsvb.b$beta, groups, ucnames)


c(
cat.res(d.test$y, as.numeric(gsvb.d.pred$mean), 
	f.gsvb.d$beta * (f.gsvb.d$g[groups] > 0.5), groups),
cat.res(d.test$y, as.numeric(gsvb.b.pred$mean), 
	f.gsvb.b$beta * (f.gsvb.b$g[groups] > 0.5), groups),
cat.res(d.test$y, as.numeric(ssgl.pred), ssgl.beta, groups)
)
