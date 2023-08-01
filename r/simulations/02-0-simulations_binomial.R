library(spsl)
library(gsvb)

source("00-functions.R")

DGP <- read.env("DGP", 1:4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 2:4)
CORES <- read.env("CORES", 1)


# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(400, 400) [SIM]
p <- c(1e3, 1e3) [SIM]
g <- c(5,   5)  [SIM]
s <- c(3,   5)  [SIM]
bmax <- 1.0
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
	list(model="binomial", corr=1),
	list(model="binomial", corr=0.6),
	list(model="binomial", corr=0.6, block_size=50),
	list(model="binomial", dof=3, weight=0.9)
    )
)

m <- list(
    # methods
    m=c(
	m_gsvb,  # GSVB (ours, jensens) 
	m_gsvb,  # GSVB (ours, jaakkola)
	m_gsvb,  # GSVB (ours, jaakola non diag)
	m_spsl   # SpSL (mcmc)
    ),
    p=list(
	list(family="binomial-jensens",  lambda=1, a0=1, b0=p/g,
	     diag_covariance=TRUE, intercept=TRUE),
	list(family="binomial-jaakkola", lambda=1, a0=1, b0=p/g,
	     diag_covariance=TRUE, intercept=TRUE),
	list(family="binomial-jaakkola", lambda=1, a0=1, b0=p/g,
	     diag_covariance=FALSE, intercept=TRUE),
	list(family="binomial", lambda=1, a0=1, b0=p/g, 
	     mcmc_samples=1e5, burnin=5e4, intercept=TRUE)
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
	     file=sprintf("../../rdata/simulations/binomial/mcmc/%s.RData", rname))
    }
}

