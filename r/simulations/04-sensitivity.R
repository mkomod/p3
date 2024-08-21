# check the sensitivity of the method wrt. to the choice of ordering in 
# the updates
library(gsvb)
source("00-functions.R")

DGP <- read.env("DGP", 1:4)
SIM <- read.env("SIM", 1)
MET <- read.env("MET", 1:6)
CORES <- read.env("CORES", 4)

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
        m_gsvb,  # GSVB (ours, with non-diagonal covariance)
        m_gsvb,  # GSVB (ours, with non-diagonal covariance)
        m_gsvb,  # GSVB (ours, with non-diagonal covariance)
        m_gsvb  # GSVB (ours, with non-diagonal covariance)
    ),
    p=list(
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=0),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=0),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=1),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=1),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=TRUE, intercept=TRUE, ordering=2),
        list(family="gaussian", lambda=1, a0=1, b0=p/g + 1, a_t=1e-3, b_t=1e-3,
             diag_covariance=FALSE, intercept=TRUE, ordering=2)
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



# n <- c(200, 200, 200, 200)
# p <- c(1e3, 1e3, 1e3, 1e3)
# g <- c(5,   5,   10,  10 )
# s <- c(5,   10,  5,   10 )
# dat <- get_data("gaussian/ordering", 
# 		1:6, 
# 		n, p, g, s,
# 		metrics=c("l2", "auc", "coverage.non_zero", "length.non_zero", "elapsed"),
# 		dgp=c(1, 2), # 1:4
# 		simnum=1, 
# 		make.table=T,
# 		make.plot=F,
# 		# fname="../figs/gaus_mcmc_1.pdf",
# 		max_psrf=3.50,
# 		method.names=c("GSVB-D", "GSVB-B", "MCMC"),
# 		method.cols=adjustcolor(color_palette[c(1,2,4)], 0.5),
# 		metric.title=c(latex2exp::TeX("$l_2$-error"), "AUC",
# 			       latex2exp::TeX("Coverage $\\beta_0 \\neq 0$"), 
# 			       latex2exp::TeX("Lenght $\\beta_0 \\neq 0$")))