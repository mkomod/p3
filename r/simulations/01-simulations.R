library(spsl)
library(gsvb)
library(sparseGAM)                      # install.packages("sparseGAM")

source("00-functions.R")

DGP <- read.env("DGP", 1:3)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 1)
CORES <- read.env("CORES", 20)


# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(200) [SIM]
p <- c(1e3) [SIM]
g <- c(5)   [SIM]
s <- c(3)   [SIM]
runs <- 100

dg <- list(
    # data generation process (dgp)
    p=c(
	dgp_diag,
	dgp_diag,
	dgp_block
    ),
    # settings for each process
    s=list(
	list(corr=0),
	list(corr=0.6),
	list(corr=0.6, block_size=50)
    )
)

m <- list(
    # methods
    m=c(
	m_gsvb,  # GSVB (ours) 
	m_gsvb,  # GSVB (ours, with non-diagonal covariance)
	m_spsl,  # SpSL (mcmc)
	m_ssgl   # SSGL (SpSL Group LASSO)
    ),
    p=list(
	list(lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3, diag_covariance=TRUE),
	list(lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3, diag_covariance=FALSE),
	list(lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3, mcmc_samples=10e3),
	list(l0=100, l1=1, a0=1, b0=p/g)
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
	rname <- sprintf("c_%d_%d_%d", i, SIM, j)
	assign(rname, m_run(m$m[[j]], m$p[[j]], setting_parameters, CORES))
	save(list=c(rname), file=sprintf("../../rdata/simulations/01/%s.RData", rname))
    }
}



