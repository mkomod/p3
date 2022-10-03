library(spsl)
library(gsvb)
library(sparseGAM)                      # install.packages("sparseGAM")

source("00-functions.R")

DGP <- read.env("DGP", 1:4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 1:4)
CORES <- read.env("CORES", 20)

if (SIM <= 4 && 3 %in% MET) MET <- MET[-which(MET == 3)]

# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(200, 200, 200, 200, 500, 500, 500, 500, 250, 1000, 250, 1000) [SIM]
p <- c(1e3, 1e3, 1e3, 1e3, 5e3, 5e3, 1e4, 1e4, 5e3,  5e3, 5e3,  5e3) [SIM]
g <- c(5,   5,   10,  10,   10,  10,  10,  10,  10,   10,  10,   10)  [SIM]
s <- c(5,   10,  5,   10,   10,  20,  10,  20,  10,   10,  20,   20)  [SIM]
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
	list(corr=0),
	list(corr=0.6),
	list(corr=0.6, block_size=50),
	list(dof=3, weight=0.9)
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
	list(family="linear", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3, 
	     mcmc_samples=10e3),
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

