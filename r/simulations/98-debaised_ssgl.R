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
n = dim(d$X)[1]
p = dim(d$X)[2]
cf = matrix(NA, p, p)
t  = rep(0, p)
S = (t(d$X) %*% d$X) / n


fit = debias_ssgl(d, fit)
ssgl_credible_intervals(fit)
method_coverage(d, fit, "ssgl")



library(gsvb)
library(sparseGAM)                      # install.packages("sparseGAM")

source("00-functions.R")

DGP <- read.env("DGP", 1:4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 1:3)
CORES <- read.env("CORES", 100)


# ----------------------------------------
# Simulation settings
# ----------------------------------------
#       1    2    3    4   5    6     7    8    9   10    11 
n <- c(250, 500, 1e3, 250, 500, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3) [SIM]
p <- c(5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3) [SIM]
g <- c( 10,  10,  10,  10,  10,  10,  10, 10,  10,  10,  10) [SIM]
s <- c( 10,  10,  10,  20,  20,  20,  10, 20,  25,  50,  100) [SIM]
bmax <- 1.5
runs <- 100

# smoke
n = 100
p = 300
g = 5
s = 5
runs = 3
CORES = 1
MET =3

dg <- list(
    # data generation process (dgp)
    p=c(
      dgp_diag,
      dgp_diag,
      dgp_block,
      dgp_wishart
    ),
    # settings for each process
    s=list(
      list(model="gaussian", corr=0),
      list(model="gaussian", corr=0.6),
      list(model="gaussian", corr=0.6, block_size=50),
      list(model="gaussian", dof=3, weight=0.9)
    )
)

m <- list(
    # methods
    m=c(
        m_gsvb,  # GSVB (ours) 
        m_gsvb,  # GSVB (ours, with non-diagonal covariance)
        m_ssgl   # SSGL (SpSL Group LASSO)
    ),
    p=list(
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3, 
             diag_covariance=TRUE, intercept=TRUE, ordering=0),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3, 
           diag_covariance=FALSE, intercept=TRUE, ordering=0),
        list(family="gaussian", l0=100, l1=1, a0=1, b0=p/g)
    )
)


# ----------------------------------------
# Run simulations
# ----------------------------------------
for (i in DGP)
{
    setting_parameters <- list(n=n, p=p, g=g, s=s, bmax=bmax, dgp=dg$p[[i]], 
	pars=dg$s[[i]], runs=runs)
    
    for (j in MET) {
        rname <- sprintf("%d_%d_%d", i, SIM, j)
        assign(rname, m_run(m$m[[j]], m$p[[j]], setting_parameters, CORES))
        # save(list=c(rname), 
	      #   file=sprintf("../../rdata/simulations/gaussian/comp/%s.RData", rname))
    }
}
