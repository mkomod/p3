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
# 			Create the validaton data
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

d.test = data.frame(y=test.Y, test.X)

# Sub-sample the data
set.seed(1234)
keep0 = sort(sample(which(d.test$y == 0), 1000, replace=FALSE))
keep1 = sort(sample(which(d.test$y == 1), 1000, replace=FALSE))

d.balanced = d.test[c(keep0, keep1), ]

# create the model matrix
X = model.matrix(y ~ .^3, d.balanced)
y = d.test$y
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




# -----------------------------------------------------------------
#   Fit GSVB
f.gsvb.d = gsvb::gsvb.fit(d$y, d$X, d$groups, family="binomial-jaakkola",
	    diag_covariance=TRUE, intercept=FALSE, track_elbo=FALSE, niter=1000)

f.gsvb.b = gsvb::gsvb.fit(d$y, d$XX, d$groups, family="binomial-jaakkola",
	    diag_covariance=FALSE, intercept=FALSE, track_elbo=FALSE, niter=1000)



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

table(d.test$y, glas.pred > 0.5)
table(d.test$y, ssgl.pred > 0.5)
table(d.test$y, gsvb.d.pred$mean > 0.5)
table(d.test$y, gsvb.b.pred$mean > 0.5)

