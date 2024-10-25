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
CORES <- read.env("CORES", 25)

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


# ----------------------------------------
# Plot
# ----------------------------------------
n <- rep(200, 15)
p <- rep(1e3, 15)
g <- rep(10, 15)
s <- rep(5, 15)
color_palette = c("#4DAF4A", "#E41A1C", "#377EB8", "#FF7F00")


METHOD = 2
dat <- get_data("gaussian/bmax", METHOD,
		n, p, g, s,
		metrics=c("l2", "auc"),
		dgp=1:4,
		simnum=1:15,
		make.table=T,
		make.plot=F,
		)

res = data.frame(); res_s = data.frame(); res_l = data.frame(); res_u = data.frame()
for (i in 1:60)
{
    # res = rbind(res, colMeans(dat[((i-1)*100 + 1):(i*100) ,]))
    res = rbind(res, apply(dat[((i-1)*100 + 1):(i*100) ,], 2, mean))
    res_s = rbind(res_s, apply(dat[((i-1)*100 + 1):(i*100) ,], 2, sd) * 1.96 / sqrt(100)) 
}

colnames(res) = colnames(dat)
colnames(res_s) = colnames(dat)
bmax <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0) 

pdf(sprintf("~/proj/papers/gsvb_ppr_ba/figures/gaus_bmax_%d.pdf", METHOD), width=12, height=5)
par(mfrow=c(1,2), mar=c(5, 5, 4, 2) + 0.1, family="Times")

m = res[res["d" ] == 1, ]$auc
s = res_s[res["d" ] == 1, ]$auc
plot(bmax, m, type="b", ylim=c(0.5, 1.0), 
    xlab=latex2exp::TeX("$\\beta_{\\max}$"), 
    ylab="AUC", 
    main="GSVB-B AUC",
    col=color_palette[1],
    pch=19, lwd=1)
polygon(c(bmax, rev(bmax)), c(m - s, rev(m + s)), col=adjustcolor(color_palette[1], 0.1), border=NA)

for (i in 2:4)
{
    m = res[res["d" ] == i, ]$auc
    s = res_s[res["d" ] == i, ]$auc
    lines(bmax, m, type="b", col=color_palette[i], pch=19, lwd=1)
    polygon(c(bmax, rev(bmax)), c(m-s, rev(m+s)), col=adjustcolor(color_palette[i], 0.1), border=NA)
}
grid()

m = res[res["d" ] == 1, ]$l2
s = res_s[res["d" ] == 1, ]$l2

plot(bmax, m, type="b", ylim=c(0, 2.75),
    xlab=latex2exp::TeX("$\\beta_{\\max}$"), 
    ylab=latex2exp::TeX("$l_2$-error"), 
    main="GSVB-B l_2-error",
    col=color_palette[1], pch=19, lwd=1)
polygon(c(bmax, rev(bmax)), c(m - s, rev(m + s)), col=adjustcolor(color_palette[1], 0.1), border=NA)

for (i in 2:4)
{
    m = res[res["d" ] == i, ]$l2
    s = res_s[res["d" ] == i, ]$l2
    lines(bmax, m, type="b", col=color_palette[i], pch=19, lwd=1)
    polygon(c(bmax, rev(bmax)), c(m-s, rev(m+s)), col=adjustcolor(color_palette[i], 0.1), border=NA)
}
grid()
dev.off()

