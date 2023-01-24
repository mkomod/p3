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

# ------------------------------------------------------------------------------
# 			 Gaussian - MCMC - GSVB
# ------------------------------------------------------------------------------
library(lattice)

n <- c(200, 200, 200, 200)
p <- c(1e3, 1e3, 1e3, 1e3)
g <- c(5,   5,   10,  10 )
s <- c(5,   10,  5,   10 )

metrics <- c("l2", "l1", "tpr", "fdr", "auc", "coverage.non_zero", 
	     "length.non_zero", "coverage.zero", "length.zero", 
	     "coverage", "elapsed")

metrics <- c("l2", "auc", "coverage.non_zero", 
	     "length.non_zero", "coverage.zero", "length.zero", 
	     "coverage", "elapsed")
dat <- c()

for (i in 1:4) { # data generating process
    for (j in 1:4) { # simulation parameters
	for (m in 1:3) { # method

	    load(file=sprintf("../../rdata/simulations/gaussian/mcmc/%d_%d_%d.RData", i, j, m))
	    x <- get(sprintf("%d_%d_%d", i, j, m))
	    if (all(c(i, j, m) == 1)) dnames <- colnames(x)
	    colnames(x) <- dnames

	    # if (sum(is.na(x[ , 1])) == 100) 
	    # {
		# next
	    # }
	    if (any(is.na(x[ , 1]))) {
		# cat("Num missing from", i, j, m, ":", sum(is.na(x[ , 1])), "\n")
	    }
	    x <- x[!is.na(x[ , 1]), ]
	    tobind <- cbind(d=i, n=n[j], p=p[j], g=g[j], s=s[j], m=m, x[ , metrics])
	    dat <- rbind(dat, tobind)

	    x.m <- t(apply(x[ , metrics], 2, function(x) 
		    quantile(x, probs=c(0.5, 0.05, 0.95), na.rm=T)))

	    for (met in seq_along(metrics)) {
		if (metrics[met] == "elapsed") {
		    times <- sapply(x.m[met, ], proc_time)
		    cat(sprintf("%s (%s, %s)", times[1], times[2], times[3]))
		} else {
		    cat(sprintf(" %.3f (%.2f, %.2f)", x.m[met, 1], x.m[met, 2], x.m[met, 3]))
		}
		cat(ifelse(met == length(metrics), "\\\\", "&"))
	    }
	    cat("\n")
	}
	cat("\n")
    }
    cat("\n\n")
}


colnames(dat)

# png("../figs/gaus.png", width=16, height=12, units="in", res=200)
{pdf("../figs/mcmc_res_gaus.pdf", width=16, height=9)
    layout(matrix(1:16, nrow=4, byrow=T))
    for (i in 1:4) 
    {
	gsize <- g[i]; n0g <- s[i]
	for (met in c("l2", "auc", "coverage.non_zero", "coverage")) 
	{
	    ddat <- dat[ dat[ , "g"] == gsize & dat[ , "s"] == n0g, ]
	    f <- as.formula(sprintf("%s ~ m*d", met))
	    # bxp
	    par(mar=c(3, 3, 1, 0))
# main=sprintf("(s=%d, g=%d)", n0g, gsize), 
	    boxplot(f, data=ddat, boxwex=0.2, col=c(2,3,4))
	}
    }
dev.off()}

metric <- c("l2", "auc", "coverage.non_zero", "coverage")
metric.title <- c("l2-error", "AUC", "Coverage non-zero", "PP Coverage")
for (i in 1:4) 
{
    pdf(sprintf("../figs/mcmc_res_gaus_%d.pdf", i), width=16, height=3)
    layout(matrix(1:4, nrow=1, byrow=T))
    gsize <- g[i]; n0g <- s[i]
    for (m in 1:4 ) 
    {
	met <- metric[m]
	ddat <- dat[ dat[ , "g"] == gsize & dat[ , "s"] == n0g, ]
	f <- as.formula(sprintf("%s ~ m*d", met))
	par(mar=c(3, 3, 3, 0))
	boxplot(f, data=ddat, main=metric.title[m], boxwex=0.2, col=c(2,3,4))
    }
    dev.off()
}




# ------------------------------------------------------------------------------
# 			 Gaussian - MCMC - GSVB
# ------------------------------------------------------------------------------
library(lattice)

n <- c(250, 500, 1e3, 250, 500, 1e3)
p <- c(5e3, 5e3, 5e3, 5e3, 5e3, 5e3)
g <- c( 10,  10,  10,  10,  10,  10)
s <- c( 10,  10,  10,  20,  20,  20)

n <- c(250, 250, 500)
p <- c(5e3, 5e3, 5e3)
g <- c( 10,  10,  10)
s <- c( 10,  20,  20)

metrics <- c("l2", "auc", "elapsed")
dat <- c()

simnum <- c(1, 4, 5)

for (j in seq_along(simnum)) { # simulation parameters
    for (i in 1:4) { # data generating process
	for (m in 1:3) { # method
	    
	    sim <- simnum[j]
	    load(file=sprintf("../../rdata/simulations/gaussian/comp/%d_%d_%d.RData", 
			      i, sim, m))
	    x <- get(sprintf("%d_%d_%d", i, sim, m))

	    if (sum(is.na(x[ , 1])) >= 90) {
		cat("NA\\\\\n")
		next
	    }

	    if (all(c(i, sim, m) == 1)) dnames1 <- colnames(x)
	    if (all(c(i, sim, m) == c(1, 1, 3))) dnames2 <- colnames(x)
	    if (any(is.null(colnames(x))) && (m == 1 || m == 2)) colnames(x) <- dnames1
	    if (any(is.null(colnames(x))) && m == 3) colnames(x) <- dnames2

	    if (any(is.na(x[ , 1]))) {
		# cat("Num missing from", i, j, m, ":", sum(is.na(x[ , 1])), "\n")
	    }
	    x <- x[!is.na(x[ , 1]), ]
	    tobind <- cbind(d=i, n=n[j], p=p[j], g=g[j], s=s[j], m=m, x[ , metrics])
	    dat <- rbind(dat, tobind)

	    x.m <- t(apply(x[ , metrics], 2, function(x) 
		    quantile(x, probs=c(0.5, 0.05, 0.95), na.rm=T)))

	    for (met in seq_along(metrics)) {
		if (metrics[met] == "elapsed") {
		    times <- sapply(x.m[met, ], proc_time)
		    cat(sprintf("%s (%s, %s)", times[1], times[2], times[3]))
		} else {
		    cat(sprintf(" %.3f (%.2f, %.2f)", x.m[met, 1], x.m[met, 2], x.m[met, 3]))
		}
		cat(ifelse(met == length(metrics), "\\\\", "&"))
	    }
	    cat("\n")
	}
	cat("\n")
    }
    cat("\n\n\n")
}


colnames(dat)

# {png("../figs/comp_gaus.png", width=12, height=12, units="in", res=200)
{pdf("../figs/comp_gaus.pdf", width=12, height=12)
    layout(matrix(1:9, nrow=3, byrow=T))
    for (i in c(1,2,3)) 
    {
	gsize <- g[i]; n0g <- s[i]; nn <- n[i]
	for (met in c("l2", "auc", "elapsed")) 
	{
	    ddat <- dat[ dat[ , "g"] == gsize & dat[ , "s"] == n0g 
			& dat[ , "n"] == nn, ]
	    f <- as.formula(sprintf("%s ~ m*d", met))

	    par(mar=c(3, 3, 1, 0))
	    boxplot(f, data=ddat, boxwex=0.2, col=c(2,3,4))
	}
    }
dev.off()}

metric <- c("l2", "auc", "elapsed")
metric.title <- c("l2-error", "AUC", "Runtime (s)")
for (i in 1:3) 
{
    pdf(sprintf("../figs/comp_gaus_%d.pdf", i), width=12, height=3)
    layout(matrix(1:3, nrow=1, byrow=T))
    gsize <- g[i]; n0g <- s[i]; nn <- n[i]
    for (m in 1:3) 
    {
	met <- metric[m]
	ddat <- dat[ dat[ , "g"] == gsize & dat[ , "s"] == n0g 
		    & dat[ , "n"] == nn, ]
	f <- as.formula(sprintf("%s ~ m*d", met))

	par(mar=c(3, 3, 3, 0))
	boxplot(f, data=ddat, main=metric.title[m], boxwex=0.2, col=c(2,3,4))
    }
    dev.off()
}


metrics <- c("coverage.non_zero", "length.non_zero", "coverage")
dat <- c()

simnum <- c(1, 4, 5)

for (j in seq_along(simnum)) { # simulation parameters
    for (i in 1:4) { # data generating process
	for (m in 1:2) { # method
	    
	    sim <- simnum[j]
	    load(file=sprintf("../../rdata/simulations/gaussian/comp/%d_%d_%d.RData", 
			      i, sim, m))
	    x <- get(sprintf("%d_%d_%d", i, sim, m))

	    if (sum(is.na(x[ , 1])) >= 90) {
		cat("NA\\\\\n")
		next
	    }

	    if (all(c(i, sim, m) == 1)) dnames1 <- colnames(x)
	    if (any(is.null(colnames(x))) && (m == 1 || m == 2)) colnames(x) <- dnames1

	    x <- x[!is.na(x[ , 1]), ]
	    tobind <- cbind(d=i, n=n[j], p=p[j], g=g[j], s=s[j], m=m, x[ , metrics])
	    dat <- rbind(dat, tobind)

	    x.m <- t(apply(x[ , metrics], 2, function(x) 
		    quantile(x, probs=c(0.5, 0.05, 0.95), na.rm=T)))

	    for (met in seq_along(metrics)) {
		if (metrics[met] == "elapsed") {
		    times <- sapply(x.m[met, ], proc_time)
		    cat(sprintf("%s (%s, %s)", times[1], times[2], times[3]))
		} else {
		    cat(sprintf(" %.3f (%.2f, %.2f)", x.m[met, 1], x.m[met, 2], x.m[met, 3]))
		}
		cat(ifelse(met == length(metrics), "\\\\", "&"))
	    }
	    cat("\n")
	}
	cat("\n")
    }
    cat("\n\n\n")
}

metric <- c("coverage.non_zero", "length.non_zero", "coverage")
metric.title <- c("Coverage Non-zero", "Length Non-zero", "PP Coverage")
for (i in 1:3) 
{
    # pdf(sprintf("../figs/cov_comp_gaus_%d.pdf", i), width=12, height=3)
    png(sprintf("../figs/cov_comp_gaus_%d.png", i), width=12, height=3, units="in", res=200)
    layout(matrix(1:3, nrow=1, byrow=T))
    gsize <- g[i]; n0g <- s[i]; nn <- n[i]
    for (m in 1:3) 
    {
	met <- metric[m]
	ddat <- dat[ dat[ , "g"] == gsize & dat[ , "s"] == n0g 
		    & dat[ , "n"] == nn, ]
	f <- as.formula(sprintf("%s ~ m*d", met))

	par(mar=c(3, 3, 3, 0))
	boxplot(f, data=ddat, main=metric.title[m], boxwex=0.2, col=c(2,3))
    }
    dev.off()
}
