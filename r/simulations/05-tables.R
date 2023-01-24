# print out the time in human read-able format
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


x <- get(sprintf("c_%d_%d_%d", 1, 1, 1))
x <- get(sprintf("%d_%d_%d", 4, 1, 1))
dnames <- colnames(x)

# ----------------------------------------
# 01-simulations
# ----------------------------------------
n <- c(200, 200, 200, 200, 250, 500, 1e3, 250, 500, 1e3, 500, 500)
p <- c(1e3, 1e3, 1e3, 1e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 1e4, 1e4)
g <- c(5,   5,   10,  10,   10,  10,  10,  10,  10,  10,  10,  10)
S <- c(5,   10,  5,   10,   10,  10,  10,  20,  20,  20,  10,  20)

all_sims <- c()

metrics <- c("l2", "l1", "tpr", "fdr", "auc", "elapsed")
for (d in 1:4) { # data generating process
    for (s in 1:11) { # simulation parameters
	for (m in 1:4) { # method

	    if (s > 4 && m == 3)
		next

	    load(file=sprintf("../../rdata/simulations/01/c_%d_%d_%d.RData", d, s, m))
	    x <- get(sprintf("c_%d_%d_%d", d, s, m))
	    if (any(m == 1:3)) colnames(x) <- dnames
	    if (sum(is.na(x[ , 1])) == 100) {
		cat("Error in run\n")
		next
	    }
	    x <- x[!is.na(x[ , 1]), ]
	    tobind <- cbind(d=d, n=n[s], p=p[s], g=g[s], s=S[s], m=m, x[ , metrics])
	    all_sims <- rbind(all_sims, tobind)

	    # x.m <- apply(x[ , metrics], 2, function(x) mean(x, na.rm=T))
	    # x.sd <- apply(x[ , metrics], 2, function(x) sd(x, na.rm=T))
	    x.m <- t(apply(x[ , metrics], 2, function(x) 
			   quantile(x, probs=c(0.5, 0.05, 0.95), na.rm=T)))

	    for (i in seq_along(metrics)) {
		if (metrics[i] == "elapsed") {
		    times <- sapply(x.m[i, ], proc_time)
		    cat(sprintf("%s (%s, %s)", times[1], times[2], times[3]))
		} else {
		    cat(sprintf(" %.3f (%.2f, %.2f)", x.m[i, 1], x.m[i, 2], x.m[i, 3]))
		}
		cat(ifelse(i == length(metrics), "\\\\", "&"))
	    }
	    cat("\n")
	}
	cat("\n")
    }
    cat("\n\n")
}


metrics <- c("coverage.non_zero", "length.non_zero", "coverage.zero", "length.zero", 
	     "coverage", "elapsed")
for (d in 1:4) { # data generating process
    for (s in 1:11) { # simulation parameters
	for (m in 1:3) { # method

	    if (s > 4 && m == 3)
		next

	    load(file=sprintf("../../rdata/simulations/01/c_%d_%d_%d.RData", d, s, m))
	    x <- get(sprintf("c_%d_%d_%d", d, s, m))
	    if (any(m == 1:3)) colnames(x) <- dnames
	    if (sum(is.na(x[ , 1])) == 100) {
		cat("Error in run\n")
		next
	    }

	    # x.m <- apply(x[ , metrics], 2, function(x) mean(x, na.rm=T))
	    # x.sd <- apply(x[ , metrics], 2, function(x) sd(x, na.rm=T))
	    x.m <- t(apply(x[ , metrics], 2, function(x) 
			   quantile(x, probs=c(0.5, 0.05, 0.95), na.rm=T)))

	    for (i in seq_along(metrics)) {
		if (metrics[i] == "elapsed") {
		    times <- sapply(x.m[i, ], proc_time)
		    cat(sprintf("%s (%s, %s)", times[1], times[2], times[3]))
		} else {
		    cat(sprintf(" %.3f (%.2f, %.2f)", x.m[i, 1], x.m[i, 2], x.m[i, 3]))
		}
		cat(ifelse(i == length(metrics), "\\\\", "&"))
	    }
	    cat("\n")
	}
	cat("\n")
    }
    cat("\n\n")
}



# ----------------------------------------
# Poisson
# ----------------------------------------
metrics <- c("l2", "l1", "tpr", "fdr", "auc", "elapsed")
for (d in 4) { # data generating process
    for (s in 1:2) { # simulation parameters
	for (m in 1:3) { # method

	    load(file=sprintf("../../rdata/simulations/poisson/%d_%d_%d.RData", d, s, m))
	    x <- get(sprintf("%d_%d_%d", d, s, m))
	    # if (any(m == 1:3)) colnames(x) <- dnames
	    if (sum(is.na(x[ , 1])) == 100) {
		cat("Error in run\n")
		next
	    }
	    x <- x[!is.na(x[ , 1]), ]
	    # tobind <- cbind(d=d, n=n[s], p=p[s], g=g[s], s=S[s], m=m, x[ , metrics])
	    # all_sims <- rbind(all_sims, tobind)

	    # x.m <- apply(x[ , metrics], 2, function(x) mean(x, na.rm=T))
	    # x.sd <- apply(x[ , metrics], 2, function(x) sd(x, na.rm=T))
	    x.m <- t(apply(x[ , metrics], 2, function(x) 
			   quantile(x, probs=c(0.5, 0.05, 0.95), na.rm=T)))

	    for (i in seq_along(metrics)) {
		if (metrics[i] == "elapsed") {
		    times <- sapply(x.m[i, ], proc_time)
		    cat(sprintf("%s (%s, %s)", times[1], times[2], times[3]))
		} else {
		    cat(sprintf(" %.3f (%.2f, %.2f)", x.m[i, 1], x.m[i, 2], x.m[i, 3]))
		}
		cat(ifelse(i == length(metrics), "\\\\", "&"))
	    }
	    cat("\n")
	}
	cat("\n")
    }
    cat("\n\n")
}



metrics <- c("l2", "l1", "tpr", "fdr", "auc",
	     "coverage.non_zero", "length.non_zero", "coverage.zero", "length.zero", 
	     "coverage", "elapsed")
for (d in 4) { # data generating process
    for (s in 1:2) { # simulation parameters
	for (m in 1:5) { # method

	    load(file=sprintf("../../rdata/simulations/binomial/%d_%d_%d.RData", d, s, m))
	    x <- get(sprintf("%d_%d_%d", d, s, m))
	    if (any(m == 1:4)) colnames(x) <- dnames
	    if (sum(is.na(x[ , 1])) == 100) {
		cat("Error in run\n")
		next
	    }
	    x <- x[!is.na(x[ , 1]), ]
	    # tobind <- cbind(d=d, n=n[s], p=p[s], g=g[s], s=S[s], m=m, x[ , metrics])
	    # all_sims <- rbind(all_sims, tobind)

	    # x.m <- apply(x[ , metrics], 2, function(x) mean(x, na.rm=T))
	    # x.sd <- apply(x[ , metrics], 2, function(x) sd(x, na.rm=T))
	    x.m <- t(apply(x[ , metrics], 2, function(x) 
			   quantile(x, probs=c(0.5, 0.05, 0.95), na.rm=T)))

	    for (i in seq_along(metrics)) {
		if (metrics[i] == "elapsed") {
		    times <- sapply(x.m[i, ], proc_time)
		    cat(sprintf("%s (%s, %s)", times[1], times[2], times[3]))
		} else {
		    cat(sprintf(" %.3f (%.2f, %.2f)", x.m[i, 1], x.m[i, 2], x.m[i, 3]))
		}
		cat(ifelse(i == length(metrics), "\\\\", "&"))
	    }
	    cat("\n")
	}
	cat("\n")
    }
    cat("\n\n")
}
