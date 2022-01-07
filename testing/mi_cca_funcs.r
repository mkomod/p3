

# compute the canonical correlation vectors for each x, y
cca <- function(x, y, objective, iter=1e3) 
{
    p <- ncol(x)
    q <- ncol(y)

    pars <- double(p + q) + rnorm(p + q)
    opt <- optim(pars, objective, x=x, y=y, control=list(fnscale=-1, maxit=iter))
    
    a <- opt$par[1:p] 
    b <- opt$par[(p+1):(q+p)]

    return(list(
	a = a / sum(sqrt(a^2)),
	b = b / sum(sqrt(b^2)),
	val = opt$value,
	iters = opt$counts[1],
	converged = !opt$convergence
    ))
}


cca.objective <- function(pars, x, y)
{
    p <- ncol(x)
    q <- ncol(y)

    a <- pars[1:p]
    b <- pars[(p+1):(q+p)]
    
    a <- a / sum(sqrt(a^2))
    b <- b / sum(sqrt(b^2))

    cor(x %*% a, y %*% b)
}


# mutual information objective
micca.objective <- function(pars, x, y)
{
    p <- ncol(x)
    q <- ncol(y)

    a <- pars[1:p]
    b <- pars[(p+1):(q+p)]

    a <- a / sum(sqrt(a^2))
    b <- b / sum(sqrt(b^2))
    
    xa <- x %*% a
    yb <- y %*% b

    px <- approxfun(density(xa))
    py <- approxfun(density(yb))

    d.xy <- MASS::kde2d(xa, yb, n=max(50, nrow(x) / 4))
    pxy <- Vectorize(
	function(x, y) d.xy$z[length(which(x >= d.xy$x)),  length(which(y >= d.xy$y))]
    )
    
    mean(log(pxy(xa, yb)) - log(px(xa)) - log(py(yb)))
}


