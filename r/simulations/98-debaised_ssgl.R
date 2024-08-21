# load in 00_functions.R
setwd("./r/simulations")
source("00-functions.R")

library(sparseGAM)
library(glmnet)

n = 100
p = 300
gsize = 5
s = 5
bmax = 1.5
pars = list(model="gaussian", corr=0)
d = dgp_diag(n, p, gsize, s, bmax, pars, b=NULL, seed=1, sig=1, n.test=100)


fit <- sparseGAM::SSGL(d$y, d$X, d$X, d$groups, family="gaussian",
lambda0=20, lambda1=1, a=1, b=100)

active_groups <- rep(0, length(unique(d$groups)))
active_groups[d$active_groups] <- 1
res <- method_summary(d$b, active_groups, fit$beta, fit$classifications, 0.5)
    
    
    
## debiasing
p0 = dim(X)[2]
cf=matrix(NA,p0,p0)
t=rep(0,p0)

sigma = (t(X) %*% X)/n
for (k in 1:p0){ 
  cvmod=cv.glmnet(x=X[,-k],y=X[,k],intercept=FALSE)
  mod=glmnet(x=X[,-k],y=X[,k],lambda=cvmod$lambda.1se,intercept=FALSE)
  cf[k,k] = 1
  cf[k,-k] = as.numeric(-mod$beta)
  
  t[k]= (1/n)*(sum((X[,k]-(X[,-k]%*%mod$beta))**2))+0.5*cvmod$lambda.1se*((sum(abs(cf[k,])))-1)
}

t1=diag(t,p0)
theta=(solve(t1))%*%cf

## debiased lasso estimate
debiased = modSSGL$beta + theta %*% t(X) %*% (Y - (X %*% modSSGL$beta))/n

## asymptotic variance of debiased lasso
varBeta =  modSSGL$sigmasq*(theta %*% sigma %*% t(theta)) / n

for (jj in 1 : (2*G)) {
  est[i,jj] = 1*(debiased[jj] - 1.96*sqrt(varBeta[jj,jj]) < beta[jj] &
                   debiased[jj] + 1.96*sqrt(varBeta[jj,jj]) > beta[jj])
  
}


