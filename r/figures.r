# ----------------------------------------
# Multivariate double exponential
# ----------------------------------------
mvde <- function(x, lambda) 
{
    x <- as.matrix(x, ncol=1)
    p <- length(x)

    x <- exp(-lambda * norm(x))
    nconst <- 2^-p * pi^(-(p-1)/2) * gamma((p+1)/2)^-1 * lambda^p

    return(nconst * x)
}


xs <- seq(-4, 4, by=0.01)
z <- outer(xs, xs, Vectorize(function(x, y) mvde(c(x, y), 1)))

nrz <- nrow(z); ncz <- ncol(z)
color <- hcl.colors(100, )
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)

jpeg("~/ap/notes/main/figures/de.jpg", width=1200, height=600)
layout(matrix(1:2, ncol=2))
par(mar=c(0, 0, 0, 0))
persp(xs, xs, z, col = color[facetcol], phi = 30, theta = -30,
      box=FALSE, border=NA, xlab="x1", ylab="x2", zlab="", useRaster=T)
par(mar=c(3, 3, 3, 3))
image(xs, xs, z, col=hcl.colors(100), useRaster=TRUE, axes=FALSE)
dev.off()

