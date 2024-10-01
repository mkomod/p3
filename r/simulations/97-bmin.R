source("./00-functions.R")
library(gsvb)

# overwrite bmin
.dgp_base <- function(n, p, gsize, s, bmax, b, seed)
{
    set.seed(seed)

    if (sum(gsize) == p) {
	groups <- rep(1:length(gsize), each=gsize)
    } else {
	if (p %% gsize != 0) stop("number of groups must be a factor of p")
	groups <- rep(1:(p/gsize), each=gsize)
    }
    
    active_groups <- sample(1:length(unique(groups)), s)

    if (is.null(b)) {
	b <- rep(0, p)
	for (a in active_groups) {
	    gk <- which(groups == a)
	    m <- length(gk)
	    b[gk] <- sample(c(1, -1), m, replace=T) * runif(m, min=0.0, max=bmax)
	}
    }
    
    return(list(groups=groups, active_groups=active_groups, b=b))
}


DGP <- read.env("DGP", 1:4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 1:2)
CORES <- read.env("CORES", 4)

# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- 200
p <- 1e3
g <- 5
s <- 10
bmax <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0) [SIM]
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
        m_gsvb
    ),
    p=list(
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=0, init_method="lasso"),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=0, init_method="lasso")
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
             file=sprintf("../../rdata/simulations/gaussian/bmax/%s.RData", rname))
    }
}

