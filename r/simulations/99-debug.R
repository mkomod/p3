# Misc code for debugging
source("./00-functions.R")

d <- dgp_diag(200, 1000, 5, 10, list(corr=0))
m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
f <- m_gsvb(d)
f <- gsvb::gsvb.fit(d$y, d$X, d$groups, intercept=TRUE, a0=1, b0=200, track_elbo=FALSE)

d <- dgp_block(200, 1000, 5, 3, list(corr=0.6, block_size=50), seed=84)
m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
f <- m_gsvb(d, m_par=m_par)

d <- dgp_block(200, 1000, 5, 3, list(corr=0.6, block_size=50), seed=91)
m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
f <- m_gsvb(d, m_par=m_par)

f <- function(x) 
{
    a <- tryCatch({
	ftime <- system.time(x)
	TRUE
    }, error=function(e) {
	return(FALSE)
    })

    if (!a) {
	return(NA)
    }
    return(ftime)
}

f(stop("123"))
f(2)
