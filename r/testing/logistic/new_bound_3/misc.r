compute_S <- function(X, m, s, g, groups) 
{
    S <- rep(1, nrow(X))

    for (group in unique(groups)) {
	G <- which(groups == group)
	S <- S * compute_S_G(X, m, s, g, G)
    }
    return(S)
}

compute_S_G <- function(X, m, s, g, G)
{
    apply(X[ , G], 1, function(x) {
	(1 - g[G][1]) + g[G][1] * exp(sum(x * m[G] + 0.5 * x^2 * s[G]^2))
    })
}

n_mgf <- function(X, m, s)
{
    apply(X, 1, function(x) {
	exp(sum(x * m + 0.5 * x^2 * s^2))
    })
}

opt_g <- function(y, X, m, s, g, G, lambda) 
{
    mk <- length(G)
    Ck <- mk * log(2) + (mk -1)/2 * log(pi) + lgamma( (mk + 1) / 2)
    S <- compute_S(X, m, s, g, groups)
    S1 <- S * n_mgf(X[ , G], m[G], s[G])

    res <- 
	log(w / (1- w)) + 
	0.5 * mk + 
	Ck +
	mk * log(lambda) +
	0.5 * sum(log(2 * pi * s[G]^2)) -
	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) +
	sum(y * X[ , G] %*% m[G]) -
	sum(log1p(S1)) +
	sum(log1p(S))

    sigmoid(res)
}


nb3 <- function(l, mu, sig)
{
    sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * pnorm(mu/sig) +
    exp(mu + 0.5*sig^2) * pnorm(-mu/sig - sig) +
    exp(-mu + 0.5*sig^2) * pnorm(mu/sig - sig)
}


xs <- seq(0, 10, by=0.05)
plot(log(1+exp(xs)), type="l")

log(1+exp(xs)) - (log(2) + xs)




