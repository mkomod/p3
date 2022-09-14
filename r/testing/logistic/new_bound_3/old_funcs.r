# ----------------------------------------
# Old functions now implemented in C++
# ----------------------------------------
e_ll <- function(X.m, X.s, ug, tau, l=20)
{
    mid <- which(ug >= tau & ug <= (1-tau))
    big <- which(ug > (1-tau))
    gk <- ug[mid]
    tot <- 0

    if (length(mid) >= 10) print("shittles")
    
    if (length(mid) == 0) {
	if (length(big) == 1) {
	    mu <- X.m[ , big]
	    sig <- sqrt(X.s[ , big])
	} else {
	    mu <- apply(X.m[ , big], 1, sum)
	    sig <- sqrt(apply(X.s[ , big], 1, sum))
	}
	return(tll(mu, sig, l))
    }

    for (i in 0:(2^length(mid)-1)) 
    {
	sk <- as.integer(intToBits(i)[1:length(mid)])
	J <- c(mid[!!sk], big)
	
	if (length(J) == 1) {
	    mu <- X.m[ , J]
	    sig <- sqrt(X.s[ , J])
	} else {
	    mu <- apply(X.m[ , J], 1, sum)
	    sig <- sqrt(apply(X.s[ , J], 1, sum))
	}

	res <- tll(mu, sig, l)
	tot <- tot + prod(gk^sk * (1 - gk)^(1 - sk)) * res
    }

    tot
}

