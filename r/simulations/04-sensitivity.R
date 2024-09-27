# check the sensitivity of the method wrt. to the choice of ordering in 
# the updates
library(gsvb)
source("00-functions.R")

DGP <- read.env("DGP", 1:4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 1:18)
CORES <- read.env("CORES", 4)

# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(200, 200, 200, 200, 200, 200, 100,  200) [SIM]
p <- c(1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 200,  200) [SIM]
g <- c(5,   5,   10,  10 , 1,   1  , 1   , 1)   [SIM]
s <- c(5,   10,  5,   10 , 5,   10 , 20  , 20)  [SIM]
bmax <- 1.5
runs <- 100

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
        m_gsvb,  
        m_gsvb,  
        m_gsvb,  
        m_gsvb,  
        m_gsvb,  
        m_gsvb,
        m_gsvb,  
        m_gsvb,  
        m_gsvb,  
        m_gsvb,  
        m_gsvb,  
        m_gsvb,
        m_gsvb,  
        m_gsvb,  
        m_gsvb,  
        m_gsvb,  
        m_gsvb,  
        m_gsvb   
    ),
    p=list(
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=0, init_method="lasso"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=0, init_method="lasso"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=1, init_method="lasso"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=1, init_method="lasso"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=2, init_method="lasso"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=2, init_method="lasso"),
        # random init   
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=0, init_method="random"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=0, init_method="random"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=1, init_method="random"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=1, init_method="random"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=2, init_method="random"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=2, init_method="random"),
       # ridge init     
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=0, init_method="ridge"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=0, init_method="ridge"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=1, init_method="ridge"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=1, init_method="ridge"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=2, init_method="ridge"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=2, init_method="ridge")
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
        save(list=c(rname), 
             file=sprintf("../../rdata/simulations/gaussian/ordering/%s.RData", rname))
    }
}

