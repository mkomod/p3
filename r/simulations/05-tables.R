# ------------------------------------------------------------------------------
library(vioplot) # install.packages("vioplot")

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
    metric.title=metrics, method.names=NULL, method.cols=methods+1,
    fname="") 
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

		if (sum(is.na(x[ , 1])) >= 95) {
		    if (make.table) cat("NA\\\\\n")
		    # x <- x[is.na(x[ , 1]), ]
		    # next
		}

		colnames(x) <- cnames[[meth]]

		if (!all(metrics %in% cnames[[meth]])) {
		    tobind = matrix(NA, nrow=nrow(x), 
			ncol=(sum(! (metrics %in% cnames[[meth]]))))
		    tobind[1, ] = 1
		    colnames(tobind) = metrics[! (metrics %in% cnames[[meth]])]
		    x = cbind(x, tobind)
		}

		tobind <- cbind(d=dtype, n=n[sim], p=p[sim], g=g[sim], 
				s=s[sim], m=meth, x[ , metrics])
		dat <- rbind(dat, tobind)

		x.m <- t(apply(x[ , metrics], 2, function(x) 
			quantile(x, probs=c(0.5, 0.05, 0.95), na.rm=T)))
		if (make.table) {
		    if (sum(is.na(x[ , 1])) >= 95) next

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
		}
		if (make.table) cat("\n")
	    }
	    if (make.table) cat("\n")
	}
	if (make.table) cat("\n\n")
    }
    
    if (make.plot) 
    {
	if (fname != "") {
	    pdf(fname, width=3.8*length(metrics), height=2.8*length(simnums))
	}
	layout(matrix(1:(length(metrics) * length(simnums)), 
		      ncol=length(metrics), byrow=T))
	for (sim in simnums)
	{
	    nn <- n[sim]; pp <- p[sim]; gg <- g[sim]; ss <- s[sim]
	    for (met in metrics)
	    {
		par(family="Times")
		if (sim == simnums[1]) {
		    par(mar=c(2, 2, 2.5, 2))
		} else {
		    par(mar=c(2, 2, 2, 2))
		}

		ddat <- dat[  dat[ , "g"] == gg & dat[ , "s"] == ss 
			    & dat[ , "n"] == nn & dat[ , "p"] == pp, ]

		f <- as.formula(sprintf("%s ~ m*d", met))

		if (met == "coverage.non_zero") {
		    rng = c(0, 1)
		} else if (met == "length.non_zero") {
		    rng = range(ddat[ , met], na.rm=T)
		    rng[1] = 0
		} else if (met == "l2") {
		    rng = range(ddat[ , met], na.rm=T)
		    rng[1] = 0
		} else if (met == "auc") {
		    # rng = range(ddat[ , met], na.rm=T)
		    rng = c(0.5, 1)
		} else {
		    rng = range(ddat[ , met], na.rm=T)
		}

		plot.new()
		plot.window(xlim=c(0.5, length(methods) * length(dgp) +0.5), 
			    ylim=rng)
		if (sim == simnums[1]) {
		    main.title = metric.title[met == metrics]
		    title(main=main.title, cex.main=2, font.main=2)
		}
		grid()
		
		# create the background alternating cols
		nmet = length(methods)
		# cols = colorRampPalette(c("white", "darkblue"))(4)[dgp]
		cols = colorRampPalette(c("white", "darkgreen"))(4)[dgp]
		# cols = c("green", "yellow", "orange", "red")[dgp]

		spacing = seq_along(dgp[-1]) * nmet + 0.5
		rect(c(-1000, spacing), -1000, c(spacing, 1000), 1000, 
		     col=adjustcolor(cols, 0.15), border=NA)
		
		# select which things to plot, if too many missing then dont plot
		# it
		att = aggregate(f, ddat, function(x) sum(is.na(x)), na.action=NULL)
		attx = 1:nrow(att)
		attx[which(att[ , 3] >= 95)] = -5
		
		# create the plot
		vioplot::vioplot(f, data=ddat, add=TRUE, wex=0.7,
		    col=method.cols, las=1, lwd=0.5, at=attx,
		    colMed=1, colMed2="white", pchMed=23, cex=1.5)

		# add the axis
		axis(1, at=1:(length(methods) * length(dgp)), 
		     labels=rep(method.names, length(dgp)), tick=T,
		    lwd=0, lwd.ticks=0.2, cex.axis=0.8)
		axis(2, at=pretty(rng), las=1, 
		    lwd=0, lwd.ticks=0.2)

		# add did not run cross
		if (any(att[ , 3] > 90)) {
		    points(which(att[ , 3] >= 90), 
			   rep(min(rng), length(which(att[ , 3] >= 90))),
			pch=4, col="red", cex=2.5)
		}
	    }
	}
	if (fname != "")
	    dev.off()
    }

    return(dat)
}


metrics <- c("l2", "l1", "tpr", "fdr", "auc", "coverage.non_zero", 
	     "length.non_zero", "coverage.zero", "length.zero", 
	     "coverage", "elapsed")


# ------------------------------------------------------------------------------
# 			 	Gaussian
# ------------------------------------------------------------------------------
n <- c(250, 500, 1e3, 250, 500, 1e3)
p <- c(5e3, 5e3, 5e3, 5e3, 5e3, 5e3)
g <- c( 10,  10,  10,  10,  10,  10)
s <- c( 10,  10,  10,  20,  20,  20)
color_palette = c("#4DAF4A", "#E41A1C", "#377EB8", "#FF7F00")

dat <- get_data("gaussian/comp", 
		1:3, 
		n, p, g, s,
		metrics=c("l2", "auc", "coverage.non_zero", "length.non_zero"), 
		dgp=1:4, 
		simnum=c(1,2,4,5), 
		make.table=FALSE,
		method.names=c("GSVB-D", "GSVB-B", "SSGL"),
		method.cols=adjustcolor(color_palette[1:3], 0.5),
		# fname="../figs/gaus_comp_1.pdf",
		metric.title=c(latex2exp::TeX("$l_2$-error"), "AUC",
			       latex2exp::TeX("Coverage $\\beta_0 \\neq 0$"), 
			       latex2exp::TeX("Lenght $\\beta_0 \\neq 0$"))
)

dat <- get_data("gaussian/comp", 
		1:3, 
		n, p, g, s,
		metrics=c("coverage.non_zero", "length.non_zero"),
		dgp=1:4, 
		simnum=c(1,2,4,5), 
		make.table=FALSE,
		method.names=c("GSVB-D", "GSVB-B", "SSGL"),
		method.cols=adjustcolor(color_palette[1:2], 0.5),
		metric.title=c(latex2exp::TeX("Coverage $\\beta_0 \\neq 0$"), 
			       latex2exp::TeX("Lenght $\\beta_0 \\neq 0$")))


n <- c(200, 200, 200, 200)
p <- c(1e3, 1e3, 1e3, 1e3)
g <- c(5,   5,   10,  10 )
s <- c(5,   10,  5,   10 )
dat <- get_data("gaussian/mcmc", 
		1:3, 
		n, p, g, s,
		metrics=c("l2", "auc", "coverage.non_zero", "length.non_zero"),
		dgp=1:4, 
		simnum=c(1,2), 
		make.table=FALSE,
		# fname="../figs/gaus_mcmc_1.pdf",
		method.names=c("GSVB-D", "GSVB-B", "MCMC"),
		method.cols=adjustcolor(color_palette[c(1,2,4)], 0.5),
		metric.title=c(latex2exp::TeX("$l_2$-error"), "AUC",
			       latex2exp::TeX("Coverage $\\beta_0 \\neq 0$"), 
			       latex2exp::TeX("Lenght $\\beta_0 \\neq 0$")))


# ------------------------------------------------------------------------------
# 			 	Binomial
# ------------------------------------------------------------------------------
n <- c(250, 500, 1e3, 250, 500, 1e3)
p <- c(5e3, 5e3, 5e3, 5e3, 5e3, 5e3) 
g <- c(  5,   5,   5,  10,  10,  10) 
s <- c(  5,   5,   5,  10,  10,  10) 

color_palette = c("darkorchid", "#4DAF4A", "#E41A1C", "#377EB8", "#FF7F00")
dat <- get_data("binomial/comp", 
		1:4, 
		n, p, g, s,
		# metrics=c("l2", "auc"), 
		metrics=c("l2", "auc", "coverage.non_zero", "length.non_zero"),
		dgp=3:4, 
		simnum=c(2,3,5,6),
		fname="../figs/binom_comp_1.pdf",
		make.table=FALSE,
		method.names=c("GSVB-D-J", "GSVB-D", "GSVB-B", "SSGL"),
		method.cols=adjustcolor(color_palette[1:4], 0.5),
		metric.title=c(latex2exp::TeX("$l_2$-error"), "AUC",
			       latex2exp::TeX("Coverage $\\beta_0 \\neq 0$"), 
			       latex2exp::TeX("Lenght $\\beta_0 \\neq 0$")))

n <- c(400, 400)
p <- c(1e3, 1e3)
g <- c(5,   5)
s <- c(3,   5)

# TODO: RUN THIS SIMULATION ASAP!!!
dat <- get_data("binomial/mcmc", 
		1:3, 
		n, p, g, s,
		metrics=c("l2", "auc", "coverage.non_zero", "length.non_zero"),
		dgp=2:4, 
		simnum=c(2,3,5,6), 
		make.table=FALSE,
		method.names=c("GSVB-D-J", "GSVB-B", "GSVB-D", "MCMC"),
		method.cols=adjustcolor(color_palette[c(1,2,3,5)], 0.5),
		metric.title=c(latex2exp::TeX("$l_2$-error"), "AUC",
		    latex2exp::TeX("Coverage $\\beta_0 \\neq 0$"),
		    latex2exp::TeX("Lenght $\\beta_0 \\neq 0$")))




# ------------------------------------------------------------------------------
# 			 	Poisson
# ------------------------------------------------------------------------------
n <- c(250, 500, 1e3, 250, 500, 1e3)
p <- c(5e3, 5e3, 5e3, 5e3, 5e3, 5e3)
g <- c(  5,   5,   5,  10,  10,  10)
s <- c(  3,   3,   3,   5,   5,   5)


color_palette = c("#4DAF4A", "#E41A1C", "#377EB8", "#FF7F00")
dat <- get_data("poisson/comp", 
		1:3, 
		n, p, g, s,
		metrics=c("l2", "auc", "coverage.non_zero", "length.non_zero"), 
		dgp=3:4, 
		simnum=c(1,2,4,5),
		make.table=FALSE,
		# fname="../figs/pois_comp_1.pdf",
		method.names=c("GSVB-D", "GSVB-D", "SSGL"),
		method.cols=adjustcolor(color_palette[1:3], 0.5),
		metric.title=c(latex2exp::TeX("$l_2$-error"), "AUC",
			       latex2exp::TeX("Coverage $\\beta_0 \\neq 0$"), 
			       latex2exp::TeX("Lenght $\\beta_0 \\neq 0$")))

n <- c(350, 350)
p <- c(1e3, 1e3)
g <- c(5,   5) 
s <- c(3,   5) 

dat <- get_data("poisson/mcmc", 
		1:3, 
		n, p, g, s,
		metrics=c("l2", "auc", "coverage.non_zero", "length.non_zero"), 
		dgp=2:4, 
		simnum=1:2,
		method.names=c("GSVB-D", "GSVB-D", "MCMC"),
		method.cols=adjustcolor(color_palette[c(1,2,4)], 0.5),
		metric.title=c(latex2exp::TeX("$l_2$-error"), "AUC",
		    latex2exp::TeX("Coverage $\\beta_0 \\neq 0$"),
		    latex2exp::TeX("Lenght $\\beta_0 \\neq 0$")))

