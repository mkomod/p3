library(spsl)
library(gsvb)

source("00-functions.R")

DGP <- read.env("DGP", 1:4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 3)
CORES <- read.env("CORES", 1)


# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(400, 400) [SIM]
p <- c(1e3, 1e3) [SIM]
g <- c(5,   5)   [SIM]
s <- c(2,   4)   [SIM]
bmax <- 0.45
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
	list(model="poisson", corr=0),
	list(model="poisson", corr=0.6),
	list(model="poisson", corr=0.6, block_size=50),
	list(model="poisson", dof=3, weight=0.9)
    )
)


m <- list(
    # methods
    m=c(
	m_gsvb,  # GSVB (ours) 
	m_gsvb,  # GSVB (ours) 
	m_spsl   # SpSL (mcmc)
    ),
    p=list(
	list(family="poisson", lambda=1, a0=1, b0=p/g, diag_covariance=TRUE, 
	     intercept=FALSE, ordering=0),
	list(family="poisson", lambda=1, a0=1, b0=p/g, diag_covariance=FALSE, 
	     intercept=FALSE, ordering=0),
	list(family="poisson", lambda=1, a0=1, b0=p/g, 
	     mcmc_samples=5e5, burnin=1e5, intercept=FALSE, kp_1=0.024, kp_2=20)
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
	save(list=c(rname), file=sprintf("../../rdata/simulations/poisson/mcmc/%s.RData", rname))
    }
}
