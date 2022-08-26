
# vars
y; X; m; s; g; groups



# ----------------------------------------
# compute_S
# ----------------------------------------
Si <- 0
i <- 1

for (j in 1:p) {
    for (k in 1:p) {
	if (j == k) {
	    ti <- g[j] * (s[j]^2 + m[j]^2)
	} else if (j != k && groups[j] == groups[k]) {
	    ti <- g[j] * m[j] * m[k]
	} else {
	    ti <- g[j] * g[k] * m[j] * m[k]
	}
	Si <- Si + X[i, j] * X[i, k] * ti 
    }
}

Si == compute_S(X, m, s, g, groups)[1]


# ----------------------------------------
# compute_S_G_K1
# ----------------------------------------
Si <- 0
i <- 1

for (j in G) {
    for (k in 1:p) {
	if (j == k) {
	    ti <- s[j]^2 + m[j]^2
	} else if (j != k && groups[j] == groups[k]) {
	    ti <- m[j] * m[k]
	} else {
	    ti <- 2 * g[k] * m[j] * m[k]
	}
	Si <- Si + X[i, j] * X[i, k] * ti 
    }
}

Si == compute_S_G_K1(X, m, s, g, G, Gc)[1]


# ----------------------------------------
# compute_S - compute_S_G
# ----------------------------------------
Si <- 0
i <- 1

for (j in Gc) {
    for (k in Gc) {
	if (j == k) {
	    ti <- g[j] * (s[j]^2 + m[j]^2)
	} else if (j != k && groups[j] == groups[k]) {
	    ti <- g[j] * m[j] * m[k]
	} else {
	    ti <- g[j] * g[k] * m[j] * m[k]
	}
	Si <- Si + X[i, j] * X[i, k] * ti 
    }
}

S <- compute_S(X, m, s, g, groups) - compute_S_G(X, m, s, g, G, Gc)
S[1] == Si

