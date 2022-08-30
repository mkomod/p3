# ----------------------------------------
# Multiplication in R
# ----------------------------------------
a <- runif(10)
B <- matrix(runif(100), nrow=10, ncol=10)

a * B == apply(B, 1, function(b) b * a)	 # rows
a * B == apply(B, 2, function(b) b * a)  # cols

# vector is multiplied by column if the same length

a <- runif(10)
B <- matrix(runif(120), nrow=12, ncol=10)

a * B == apply(B, 1, function(b) b * a)	 # rows
a * B == apply(B, 2, function(b) b * a)  # cols

# vector same lenght as row will not raise error if multiplied



# ----------------------------------------
# Testing compute_S
# ----------------------------------------
S <- c()
for (i in 1:n) {
    temp <- 1
    for (group in unique(groups)) {
	G <- which(groups == group)
	temp <- temp * (
	    (1 - g[G][1]) + 
	    g[G][1] * exp(sum(l[i] * X[i, G] * m[G] + 0.5 * 
			      l[i]^2 * X[i, G]^2 * s[G]^2))
	)
    }
    S <- c(S, temp)
}

# ----------------------------------------
# Testing the bound
# ----------------------------------------
# - E log s(-x) <= 
#	E lx + log( E e^(-lx) + E e^(x(1-l)) )
# s(x) = (1 + e^(-x))^-1

s <- function(x) 1/(1 +  exp(-x))

nmgf <- function(m, s, l) exp(l * m + 0.5 * l^2 * s^2)

x <- rnorm(1000)

mcs <- mean(-log(s(x)))
jen <- log(1 + nmgf(0, 1, 1))

optim(1, function(l) l * 0 + log(nmgf(0, 1, -l) + nmgf(0, 1, 1-l)))

lls <- seq(-0.2, 1.2, length.out=200)
ll <- function(l) l * 0 + log(nmgf(0, 1, -l) + nmgf(0, 1, 1-l))

plot(lls, sapply(lls, ll), type="l", ylim=c(0.5, 1))
abline(h=jen)
abline(h=mcs)

