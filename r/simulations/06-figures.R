# ------------------------------------------------------------------------------
proc_time <- function(e) 
{
    if (e < 60) {
	sprintf("%.1fs", e)
    } else if (e <= 3600) {
	m <- e / 60
	s <- (m - floor(m)) * 60
	sprintf("%.0fm %.0fs", m, s)
    } else {
	h <- e / 3600
	m <- (h - floor(h)) * 60
	sprintf("%.0fh %.0fm", floor(h), m)
    }
}


get_data <- function(family, methods, n, p, g, s, metrics,
    dgp=1:4, simnums=1:length(n), make.table=FALSE, make.plot=TRUE,
    metric.title=NULL) 
{
    cnames <- list() 
    for (i in methods) {
	load(file=sprintf("../../rdata/simulations/%s/%d_%d_%d.RData", 
			  family, 1, 1, i))
	dnames <- colnames(get(sprintf("1_1_%d", i)))
	if (is.null(dnames)) {
	    cnames[[i]] <- cnames[[i-1]]
	} else {
	    cnames[[i]] <- dnames
	}
    }

    dat <- c()

    for (dtype in dgp) { 
	for (sim in simnums) { # simulation parameters
	    for (meth in methods) 
	    {
		load(file=sprintf("../../rdata/simulations/%s/%d_%d_%d.RData", 
				  family, dtype, sim, meth))
		x <- get(sprintf("%d_%d_%d", dtype, sim, meth))

		if (sum(is.na(x[ , 1])) >= 90) {
		    if (make.table) cat("NA\\\\\n")
		    next
		}

		colnames(x) <- cnames[[meth]]

		x <- x[!is.na(x[ , 1]), ]
		tobind <- cbind(d=dtype, n=n[sim], p=p[sim], g=g[sim], 
				s=s[sim], m=meth, x[ , metrics])
		dat <- rbind(dat, tobind)

		x.m <- t(apply(x[ , metrics], 2, function(x) 
			quantile(x, probs=c(0.5, 0.05, 0.95), na.rm=T)))
		if (make.table)
		    for (met in seq_along(metrics)) {
			if (metrics[met] == "elapsed") {
			    times <- sapply(x.m[met, ], proc_time)
			    cat(sprintf("%s (%s, %s)", times[1], times[2], times[3]))
			} else {
			    cat(sprintf(" %.3f (%.2f, %.2f)", x.m[met, 1], 
					x.m[met, 2], x.m[met, 3]))
			}
			cat(ifelse(met == length(metrics), "\\\\", "&"))
		    }
		if (make.table) cat("\n")
	    }
	    if (make.table) cat("\n")
	}
	if (make.table) cat("\n\n")
    }
    
    if (make.plot) 
    {
	layout(matrix(1:(length(metrics) * length(simnums)), 
		      ncol=length(metrics), byrow=T))
	for (sim in simnums)
	{
	    nn <- n[sim]; pp <- p[sim]; gg <- g[sim]; ss <- s[sim]
	    for (met in metrics)
	    {
		mm <- metrics[met]
		ddat <- dat[  dat[ , "g"] == gg & dat[ , "s"] == ss 
			    & dat[ , "n"] == nn & dat[ , "p"] == pp, ]

		# mat <- c()
		# for (met in metrics)
		#     mat <- cbind(mat, ddat[ddat[ , "m"] == meth, met])

		# par(mar=c(3, 3, 1, 0))
		# boxplot(mat, data=ddat, boxwex=0.2, col=methods+1, drop=FALSE)
		f <- as.formula(sprintf("%s ~ m*d", met))

		par(mar=c(3, 3, 3, 0))
		boxplot(f, data=ddat, boxwex=0.2, col=methods+1)
	    }
	}
    }

    return(dat)
}


# ------------------------------------------------------------------------------
# 			 	Gaussian
# ------------------------------------------------------------------------------
library(lattice)

metrics <- c("l2", "l1", "tpr", "fdr", "auc", "coverage.non_zero", 
	     "length.non_zero", "coverage.zero", "length.zero", 
	     "coverage", "elapsed")

n <- c(250, 500, 1e3, 250, 500, 1e3)
p <- c(5e3, 5e3, 5e3, 5e3, 5e3, 5e3)
g <- c( 10,  10,  10,  10,  10,  10)
s <- c( 10,  10,  10,  20,  20,  20)
dat <- get_data("gaussian/comp", 1:3, n, p, g, s, c("l2", "auc"), dgp=3:4, simnum=c(1,2,4,5))
dat <- get_data("gaussian/comp", 1:3, n, p, g, s, c("l2", "auc"), simnum=c(4,5), make.table=T)
dat <- get_data("gaussian/comp", 1:2, n, p, g, s, c("coverage.non_zero", "length.non_zero", 
						    "coverage"), simnum=c(1,2,4,5))
dat <- get_data("gaussian/comp", 1:3, n, p, g, s, c("l2", "auc"), 
		dgp=3:4, simnum=c(1, 2), make.table=T)

n <- c(200, 200, 200, 200)
p <- c(1e3, 1e3, 1e3, 1e3)
g <- c(5,   5,   10,  10 )
s <- c(5,   10,  5,   10 )
dat <- get_data("gaussian/mcmc", 1:3, n, p, g, s, c("l2", "auc", "elapsed"))
dat <- get_data("gaussian/mcmc", 1:3, n, p, g, s, c("l2", "auc", "coverage.non_zero", 
						    "coverage"))
dat <- get_data("gaussian/mcmc", 1:3, n, p, g, s, c("l2", "auc", "coverage.non_zero", "coverage", "elapsed"), 
		dgp=3:4, simnum=4, make.table=T)

# ------------------------------------------------------------------------------
# 			 	Binomial
# ------------------------------------------------------------------------------
n <- c(250, 500, 1e3, 250, 500, 1e3)
p <- c(5e3, 5e3, 5e3, 5e3, 5e3, 5e3) 
g <- c(  5,   5,   5,  10,  10,  10) 
s <- c(  5,   5,   5,  10,  10,  10) 

dat <- get_data("binomial/comp", 1:4, n, p, g, s, c("l2", "auc", "elapsed"), dgp=2:4)
dat <- get_data("binomial/comp", 1:3, n, p, g, s, c("coverage.non_zero", "length.non_zero"), 
		dgp=4)

n <- c(350, 350)
p <- c(1e3, 1e3)
g <- c(5,   5)
s <- c(3,   5)

dat <- get_data("binomial/mcmc", 1:2, n, p, g, s, c("l2", "auc", "elapsed"), dgp=1:2)
dat <- get_data("binomial/comp", 1:3, n, p, g, s, c("coverage.non_zero", "length.non_zero"), 
		dgp=4)


# ------------------------------------------------------------------------------
# 			 	Poisson
# ------------------------------------------------------------------------------
n <- c(250, 500, 1e3, 250, 500, 1e3)
p <- c(5e3, 5e3, 5e3, 5e3, 5e3, 5e3)
g <- c(  5,   5,   5,  10,  10,  10)
s <- c(  3,   3,   3,   5,   5,   5)

dat <- get_data("poisson/comp", 1:3, n, p, g, s, c("l2", "auc", "elapsed"), dgp=1:4, simnum=1:5)
dat <- get_data("poisson/comp", 1:2, n, p, g, s, c("coverage.non_zero", "length.non_zero"), 
		dgp=4, simnum=1:5)


n <- c(350, 350)
p <- c(1e3, 1e3)
g <- c(5,   5) 
s <- c(3,   5) 
dat <- get_data("poisson/mcmc", 1:3, n, p, g, s, c("l2", "auc", "elapsed"), dgp=1:4)
dat <- get_data("poisson/mcmc", 1:3, n, p, g, s, c("coverage.non_zero", "length.non_zero"), 
		dgp=1:4)
