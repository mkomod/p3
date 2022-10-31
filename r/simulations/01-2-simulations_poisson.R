library(spsl)
library(gsvb)
library(sparseGAM)                      # install.packages("sparseGAM")

source("00-functions.R")

DGP <- read.env("DGP", 4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 1:3)
CORES <- read.env("CORES", 3)

# if (SIM > 4 && 3 %in% MET) MET <- MET[-which(MET == 3)]

# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(250, 250, 250, 250) [SIM]
p <- c(1e3, 1e3, 1e3, 1e3) [SIM]
g <- c(5,   5,   10,  10)  [SIM]
s <- c(5,   10,  5,   10)  [SIM]
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
	m_spsl,  # SpSL (mcmc)
	m_ssgl   # SSGL (SpSL Group LASSO)
    ),
    p=list(
	list(family="poisson", lambda=1, a0=1, b0=p/g + 1, diag_covariance=TRUE, 
	     intercept=FALSE),
	list(family="poisson", lambda=1, a0=1, b0=p/g + 1, mcmc_samples=10e3,
	     intercept=FALSE),
	list(family="poisson", l0=100, l1=1, a0=1, b0=p/g)
    )
)


# ----------------------------------------
# Run simulations
# ----------------------------------------
for (i in DGP)
{
    setting_parameters <- list(n=n, p=p, g=g, s=s, dgp=dg$p[[i]], 
	pars=dg$s[[i]], runs=runs)
    
    for (j in MET) {
	rname <- sprintf("%d_%d_%d", i, SIM, j)
	assign(rname, m_run(m$m[[j]], m$p[[j]], setting_parameters, CORES))
	save(list=c(rname), file=sprintf("../../rdata/simulations/poisson/%s.RData", rname))
    }
}

