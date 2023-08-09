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
cv.glasso = gglasso::cv.gglasso(d$XX, d$yy, d$groups, loss="logit",
	    pred.loss="misclass", intercept=FALSE)
plot(cv.glasso)
lambda.min = cv.glasso$lambda[which.min(cv.glasso$cvlo)] # lambda.min =  0.001766377

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

f.ssgl = sparseGAM::SSGL(d$y, d$XX, d$XX[1:10, ], groups, family="binomial",
		lambda0=cv.ssgl$lambda0.min)

# f.ssgl = sparseGAM::SSGL(y, X, X[1:10, ], groups, family="binomial",
# 		lambda0=25)



# -----------------------------------------------------------------
#   Fit GSVB
f.gsvb.d = gsvb::gsvb.fit(d$y, d$XX, d$groups, family="binomial-jaakkola", lambda=1,
	    diag_covariance=TRUE, intercept=FALSE, track_elbo=FALSE, niter=1000)

f.gsvb.b = gsvb::gsvb.fit(d$y, d$XX, d$groups, family="binomial-jaakkola", lambda=1,
	    diag_covariance=FALSE, intercept=FALSE, track_elbo=FALSE, niter=1000)

n.d = plot.norms(f.gsvb.d$beta, groups)
n.b = plot.norms(f.gsvb.b$beta, groups)
n.g = plot.norms(f.gg$beta, groups)

plot(n.d, ylim=c(0, 10))
points(n.b, pch=2)
points(n.g, pch=3)

norms.d = comp.norms(f.gsvb.d, groups)
q.norms.d = apply(norms.d, 1, quantile, probs=c(0.025, 0.975))





# -----------------------------------------------------------------------------
# 			Plot the norms of the methods
# -----------------------------------------------------------------------------
pdf("../figs/memset_norms.pdf", 9, 6)
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,3), ncol=1))

# GSVB-D
par(mar=c(0, 4, 0.1, 1), family="Times")
plot.norms(f.gsvb.d$beta, groups, ylim=c(0, 3.5), pch=20, xaxt="n", las=1,
	   col=ifelse(f.gsvb.d$g > 0.5, 1, "lightgrey"), ylab="GSVB-D")
arrows((1:max(groups)) * (-1)^((f.gsvb.d$g > 0.5) - 1), 
       q.norms.d[1, ], y1=q.norms.d[2, ], code=0, col=1, lwd=2)
selected = (f.gsvb.d$g > 0.5 & f.gsvb.b$g > 0.5) | (f.gsvb.d$g > 0.5 & n.g)
abline(v = which(selected), lty=2, lwd=0.5, col=2)
grid(nx=0, ny=5)

# GSVB-B
par(mar=c(0, 4, 0, 1))
plot.norms(f.gsvb.b$beta, groups, ylim=c(0, 13.5), pch=20, xaxt="n", las=1,
	   col=ifelse(f.gsvb.b$g > 0.5, 1, "lightgrey"), ylab="GSVB-B")
arrows((1:max(groups)) * (-1)^((f.gsvb.b$g > 0.5) - 1), 
       q.norms.b[1, ], y1=q.norms.b[2, ], code=0, col=1, lwd=2)
selected = (f.gsvb.b$g > 0.5 & f.gsvb.d$g > 0.5) | (f.gsvb.b$g > 0.5 & n.g)
abline(v = which(selected), lty=2, lwd=0.5, col=2)
grid(nx=0, ny=5)
par(mar=c(3, 4, 0, 1))


# Group LASSO
plot.norms(f.gg$beta, groups, ylim=c(0, 10.5), pch=20, las=1, xaxt="n",
	   col=ifelse(n.g, 1, "lightgrey"), ylab="Group LASSO")
selected = (n.g & f.gsvb.d$g > 0.5) | (f.gsvb.b$g > 0.5 & n.g)
abline(v = which(selected), lty=2, lwd=0.5, col=2)
axis(1, 1:max(groups), F, lwd=0, lwd.ticks=0.8)
text(x=1:max(groups), par("usr")[3] - 0.8,
     labels = ucnames, srt = -45, 
     pos = 1, xpd = TRUE, cex=0.82)
grid(nx=0, ny=5)


dev.off()








# -----------------------------------------------------------------
#   Save everything
save(list=c("f.gsvb.b", "f.gsvb.d", "f.gg", "cv.glasso", "cv.ssgl", "f.ssgl"),
     file="../../rdata/application/memset/models.RData")


# -----------------------------------------------------------------
#   Load
load("../../rdata/application/memset/models.RData")



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
glas.pred = 1 / (1 + exp(- d.test$XX %*% f.gg$beta))
# ssgl.pred = 1 / (1 + exp(- cbind(1, d.test$X) %*% c(f.ssgl$beta0[1], f.ssgl$beta)))
gsvb.d.pred = gsvb::gsvb.predict(f.gsvb.d, d.test$X)
gsvb.b.pred = gsvb::gsvb.predict(f.gsvb.b, d.test$XX)

table(d.test$y, glas.pred > 0.5)
# table(d.test$y, ssgl.pred > 0.5)
table(d.test$y, gsvb.d.pred$mean > 0.5)
table(d.test$y, gsvb.b.pred$mean > 0.5)


# compute rho_max

taus = seq(0.01, 0.99, by=0.01)
max(sapply(taus, function(tau) cor(d.test$y, glas.pred > tau)))
max(sapply(taus, function(tau) cor(d.test$y, gsvb.d.pred$mean > tau))) 
max(sapply(taus, function(tau) cor(d.test$y, gsvb.b.pred$mean > tau)))



# ROSE::roc.curve(d.test$y, as.numeric(glas.pred), col=3, add=T)
ROSE::roc.curve(d.test$y, as.numeric(glas.pred))
ROSE::roc.curve(d.test$y, as.numeric(ssgl.pred))
ROSE::roc.curve(d.test$y, as.numeric(gsvb.d.pred$mean), col=2, add=T)
ROSE::roc.curve(d.test$y, as.numeric(gsvb.b.pred$mean), col=3, add=T)

ROSE::accuracy.meas(d.test$y, as.numeric(glas.pred))
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

# correlation between true and predicted for different threshold



