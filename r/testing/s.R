xs <- seq(-50, 50, 0.1)
fx <- log(1+exp(xs))
gx <- sapply(xs, function(x) 1 + max(0, x))
hx <- 1 + 0.5 * (xs + abs(xs))
ix <- 0.5 * (xs + sqrt(2 + xs^2))

plot(xs, fx, type="l", lwd=4, col=2, xlab="x", ylab="f(x)")
lines(xs, ix, lwd=4, col=3)
plot(xs, fx - ix)
legend("bottomright", legend=c("log(1 + exp(x))", "x + sqrt(2 + x^2)"),
       col=c(2,3), lwd=4)
grid()
savePlot("~/p1/notes/gsvb/figures/bounds.jpg", type="jpeg")

xs <- rnorm(500)
mean(fx)
mean(ix)

mean(0.5 * (xs + sqrt(2 + mean(xs^2))))

xs <- seq(0, 1, 0.001)
plot(xs, xs, type="l")
lines(xs, sqrt(xs))


g <- seq(0, 1, 0.01)
fx <- sapply(g, function(g) {
    xs <- rnorm(9e4) * (runif(9e4) <= g)
    mean(sqrt(2 + xs^2))
})
plot(g, g * sqrt(2 + 5) + (1-g) * sqrt(2), type="l")
lines(g, sqrt(2 + g * 5))
lines(c(g[1], g[101]), c(fx[1], fx[101]))

xs <- seq(0, 5, 0.1)
plot(xs, sqrt(2 + xs))
plot(xs, sqrt(xs))


ys <- rnorm(500000, 2, 1) * (runif(500000) > 0.2)
xs <- rnorm(500000, -2, 1) * (runif(500000) > 0.9)

mean(sqrt(2 + xs^2 + 2*ys*xs + ys^2))

xs0 <- rep(0, 500000)
xs1 <- rnorm(500000, -2, 1)

0.1 * mean(sqrt(2 + xs1^2 + 2*ys*xs1 + ys^2)) + 
0.9 * mean(sqrt(2 + xs0^2 + 2*ys*xs0 + ys^2))

xs0 <- rep(0, 500)
xs1 <- rnorm(500, -2, 1)

0.1 * sqrt(2 + mean(xs1^2 + 2*ys*xs1 + ys^2)) + 0.9 * sqrt(2 + mean(xs0^2 + 2*ys*xs0 + ys^2))



