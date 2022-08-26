compute_S <- function(X, m, s, g, groups) 
{
    gg <- g[!duplicated(groups)]
    tmp <- matrix(1, nrow=length(gg), ncol=length(gg))
    diag(tmp) <- ifelse(gg, gg, 1)
    G <- outer(gg, gg) / tmp

    E_bi_bj <- (outer(m , m) + diag(s^2)) * G[groups, groups]

    apply(X, 1, function(x) {
	sum(outer(x, x) * E_bi_bj)
    })
}

compute_S_G <- function(X, m, s, g, G, Gc)
{
    A <- outer(m[G], m[Gc]) * outer(g[G], g[Gc])
    B <- (outer(m[G], m[G]) + diag(s[G]^2)) * g[G][1] 
    apply(X, 1, function(x) {
	2 * sum(outer(x[G], x[Gc]) * A) + sum(outer(x[G], x[G]) * B)
    })
}


compute_S_G_K1 <- function(X, m, s, g, G, Gc)
{
    A <- outer(m[G], m[Gc]) * outer(rep(1, length(G)), g[Gc])
    B <- (outer(m[G], m[G]) + diag(s[G]^2))
    apply(X, 1, function(x) {
	2 * sum(outer(x[G], x[Gc]) * A) + sum(outer(x[G], x[G]) * B)
    })
}



opt_mu <- function(m_G, y, X, m, s, g, G, Gc, lambda, S) 
{
    m[G] <- m_G
    S <- S + compute_S_G_K1(X, m, s, g, G, Gc)

    sum((0.5 - y) * (X[ , G] %*% m_G) + 0.5 * (2 + S)^0.5) -
    lambda * sqrt(sum(s[G]^2) + sum(m_G^2))
}


opt_s <- function(s_G, y, X, m, s, g, G, Gc, lambda, S) 
{
    s[G] <- s_G
    S <- S + compute_S_G_K1(X, m, s, g, G, Gc)

    sum(0.5 * sqrt(2 + S)) -
    sum(log(s[G])) +
    lambda * sqrt(sum(s[G]^2) + sum(m[G]^2))
}

opt_g <- function(y, X, m, s, g, G, Gc, lambda, S) 
{
    mk <- length(G)
    Ck <- mk * log(2) + (mk -1) * log(pi) + lgamma( (mk + 1) / 2)
    S1 <- S + compute_S_G_K1(X, m, s, g, G, Gc)

    res <- log(w / (1- w)) + 0.5 * mk + 
	0.5 * sum(log(2 * pi * s[G]^2)) -
	Ck -
	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) +
	mk * log(lambda) -
	sum(
	    (0.5 - y) * X[ , G] %*% m[G] +
	    0.5 * sqrt(2 + S1) -
	    0.5 * sqrt(2 + S)
	)
    sigmoid(res)
}


sigmoid <- function(x) 1.0 / (1.0 + exp(-x))

