library(gsvb)
library(sparseGAM)                      # install.packages("sparseGAM")

source("00-functions.R")

DGP <- read.env("DGP", 1:4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 1:2)
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
	     diag_covariance=TRUE, intercept=TRUE, ordering=2),
	list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3, 
	     diag_covariance=FALSE, intercept=TRUE, ordering=2),
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
	save(list=c(rname), 
	     file=sprintf("../../rdata/simulations/gaussian/comp/%s.RData", rname))
    }
}

