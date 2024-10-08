library(spsl)
library(gsvb)
library(sparseGAM)                      # install.packages("sparseGAM")

source("00-functions.R")

DGP <- read.env("DGP", 1:4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 4)
CORES <- read.env("CORES", 1)

# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(200, 200, 200, 200, 400, 400, 400, 400) [SIM]
p <- c(1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3) [SIM]
g <- c(5,   5,   10,  10 , 5,   5,   10,  10 ) [SIM]
s <- c(5,   10,  5,   10 , 5,   10,  5,   10 ) [SIM]
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
	m_gsvb,  # GSVB (ours) 
	m_gsvb,  # GSVB (ours, with non-diagonal covariance)
	m_spsl,  # SpSL (mcmc)
	m_ssgl 	 # SpSL (mcmc)
    ),
    p=list(
	list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
	     diag_covariance=TRUE, intercept=TRUE, ordering=0),
	list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
	     diag_covariance=FALSE, intercept=TRUE, ordering=0),
	list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
	     mcmc_samples=1e5, burnin=5e4, intercept=TRUE),
	list(family="gaussian", l0=20, l1=1, a0=1, b0=p/g)
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
	     file=sprintf("../../rdata/simulations/gaussian/mcmc/%s.RData", rname))
    }
}

