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


# ----------------------------------------
# 01-simulations
# ----------------------------------------
metrics <- c("l2", "l1", "tpr", "fdr", "auc", "elapsed")
for (d in 1:3) { # data generating process
    for (s in 1) { # simulation parameters
	for (m in 1:4) { # method

	    load(file=sprintf("../../rdata/simulations/01/c_%d_%d_%d.RData", d, s, m))
	    x <- get(sprintf("c_%d_%d_%d", d, s, m))
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
    cat("\n")
}


