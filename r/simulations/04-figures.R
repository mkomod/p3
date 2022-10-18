
x <- get(sprintf("c_%d_%d_%d", 1, 1, 1))
dnames <- colnames(x)

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
		next
	    }
	    x <- x[!is.na(x[ , 1]), ]
	    tobind <- cbind(d=d, n=n[s], p=p[s], g=g[s], s=S[s], m=m, x[ , metrics])
	    all_sims <- rbind(all_sims, tobind)
	}
    }
}

library(lattice)
all_sims <- data.frame(all_sims)
all_sims$d.f <- as.factor(all_sims$d)
all_sims$n.f <- as.factor(all_sims$n)
all_sims$p.f <- as.factor(all_sims$p)
all_sims$g.f <- as.factor(all_sims$g)
all_sims$s.f <- as.factor(all_sims$s)
all_sims$m.f <- as.factor(all_sims$m)

plt <- bwplot(l2 ~ m.f|g.f*s.f*d.f, data=all_sims[all_sims$n == 200, ])
plt[ , , 1]
plt[ , , 2]
plt[ , , 3]
plt[ , , 4]
plt <- bwplot(auc ~ m.f|g.f*s.f*d.f, data=all_sims[all_sims$n == 200, ])
plt[ , , 1]
plt[ , , 2]
plt[ , , 3]
plt[ , , 4]

plt <- bwplot(l2 ~ m.f|n.f*s.f*d.f, data=all_sims[all_sims$p == 5000, ])
plt[ , , 1]
plt[ , , 2]
plt[ , , 3]
plt[ , , 4]
plt <- bwplot(auc ~ m.f|n.f*s.f*d.f, data=all_sims[all_sims$p == 5000, ])
plt[ , , 1]
plt[ , , 2]
plt[ , , 3]
plt[ , , 4]

n <- c(200, 200, 200, 200, 250, 500, 1e3, 250, 500, 1e3, 500, 500)
p <- c(1e3, 1e3, 1e3, 1e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 1e4, 1e4)
g <- c(5,   5,   10,  10,   10,  10,  10,  10,  10,  10,  10,  10)
S <- c(5,   10,  5,   10,   10,  10,  10,  20,  20,  20,  10,  20)

all_sims <- c()

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
		next
	    }
	    x <- x[!is.na(x[ , 1]), ]
	    tobind <- cbind(d=d, n=n[s], p=p[s], g=g[s], s=S[s], m=m, x[ , metrics])
	    all_sims <- rbind(all_sims, tobind)
	}
    }
}

all_sims <- data.frame(all_sims)
all_sims$d.f <- as.factor(all_sims$d)
all_sims$n.f <- as.factor(all_sims$n)
all_sims$p.f <- as.factor(all_sims$p)
all_sims$g.f <- as.factor(all_sims$g)
all_sims$s.f <- as.factor(all_sims$s)
all_sims$m.f <- as.factor(all_sims$m)


plt <- bwplot(coverage.non_zero ~ m.f|g.f*s.f*d.f, data=all_sims[all_sims$n == 200, ])
plt[ , , 1]
plt[ , , 2]
plt[ , , 3]
plt[ , , 4]
plt <- bwplot(length.non_zero ~ m.f|g.f*s.f*d.f, data=all_sims[all_sims$n == 200, ])
plt[ , , 1]
plt[ , , 2]
plt[ , , 3]
plt[ , , 4]
plt <- bwplot(coverage ~ m.f|g.f*s.f*d.f, data=all_sims[all_sims$n == 200, ], ylim=c(0.8, 1.1))
plt[ , , 1]
plt[ , , 2]
plt[ , , 3]
plt[ , , 4]

