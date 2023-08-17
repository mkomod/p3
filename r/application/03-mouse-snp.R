# -------------------------------------------------------------------------------
# 			      Mouse Genotype Data
#
# -------------------------------------------------------------------------------
rm(list=ls())

source("./00-functions.R")


get_genes = function(snps) 
{
    keep = snps %in% snp_locations$marker
    snps = snps[keep]
    snp_locations = snp_locations[ snp_locations$marker %in% snps, ]

    marker_locations = marker_locations[marker_locations$Marker.Type == "Gene", ]

    genes = apply(snp_locations, 1, function(snp) 
    {
	# make sure the chromosome is the same
	marker_locations_sub = marker_locations[marker_locations$Chrom == snp[3], ]

	# is the location of the SNP wiithin the gene location
	bp = as.numeric(snp[4])
	# marker_locations_sub$Marker.Symbol[
	#     which(marker_locations_sub$Start < bp & bp < marker_locations_sub$End)
	# ]
	marker_locations_sub[
	    which(marker_locations_sub$Start < bp & bp < marker_locations_sub$End),
	]
    })

    genes
}


pred.y = function(fit, X.star, y.train, X.train, sample=1500)
{
    samp = gsvb::gsvb.sample(fit, sample=sample)
    inter = sweep(X.train %*% samp$beta, 1, matrix(y.train, ncol=1))
    inter = colMeans(inter)

    Xb = sweep(X.star %*% samp$beta, 2, inter, "-")

    sigma <- sqrt(fit$tau_b / fit$tau_a)
    y.star <- Xb + sigma * rt(prod(dim(Xb)), 2 + fit$tau_a)
    
    return(y.star)
}


mse = function(x, y) mean( (y - x)^2)

# -------------------------------------------------------------------------------
# 			        Load the data
# -------------------------------------------------------------------------------
genotypes = read.csv("../../data/mouse/X.csv")
phenotypes = read.csv("../../data/mouse/phenotypes.csv")
snp_locations = read.csv("../../data/mouse/mice_map2.csv")
marker_locations = read.csv("../../data/mouse/MGI_MRK_Coord.rpt.csv")


# -------------------------------------------------------------------------------
# Select the phenotype of interest and re-order the design
# -------------------------------------------------------------------------------
phenotypes = phenotypes[order(phenotypes[ , 1]), ]
genotypes = genotypes[order(genotypes[ , 1]), ]

all(genotypes[ , 1] == phenotypes[ , 1])
ids = genotypes[ , 1]

# make the 
rownames(genotypes) = ids
rownames(phenotypes) = ids
genotypes = genotypes[ , -1]
phenotypes = phenotypes[ , -1]


# Create the response
y = phenotypes$Biochem.LDL
keep = !is.na(y)

# keep non missing values
y = y[keep]
genotypes = genotypes[keep, ]


length(y)
dim(genotypes)


# -------------------------------------------------------------------------------
# 			 Pre-process the genotypes
# -------------------------------------------------------------------------------
any(is.na(genotypes))
n = length(y)

# compute the minor allel frequency
mafX = apply(genotypes, 2, function(x) {
    (sum(x == 1) + 2 * sum(x == 2)) / (2 * n)
})
mafXq25 = quantile(mafX, prob=0.25)

keep = which(mafX > mafXq25)
genotypes = genotypes[ , keep]

# varX = apply(genotypes, 2, var)
# varXq25 = quantile(varX, prob=0.25)
# keep = which(varX > varXq25)
# genotypes = genotypes[ , keep]


# clean up the genotypes
snps = gsub("_[A|T|C|G]$", "", colnames(genotypes))
colnames(genotypes) = snps

# remove the geneotypes that are in high colinearity with on another
X = as.matrix(genotypes)
S = cor(X)

p = ncol(X)
keep = rep(T, p)
in_cor = list()

for (i in 1:(p - 1)) 
{
    if (keep[i]) {
	high_cor = i + which(abs(S[i, (i + 1):p]) > 0.97)
	keep[high_cor] = F
	in_cor[[snps[i]]] = snps[high_cor]
    }
}

genotypes = genotypes[ , keep]
snps = snps[keep]

dim(genotypes)

# -------------------------------------------------------------------------------
# 			Create the model matrix
# -------------------------------------------------------------------------------
X = apply(as.matrix(genotypes), 2, as.factor)
X = data.frame(X, stringsAsFactors=T)

con2 = matrix(c(0, 1), nrow=2)
con3 = matrix(c(0, 0, 0, 1, 1, 1), nrow=3, byrow=T)

for (i in names(X)) {
    if (length(levels(X[[i]])) == 3) {
	contrasts(X[[i]]) = con3
    } else {
	contrasts(X[[i]]) = con2
    }
}

X = model.matrix(y ~ . - 1, data.frame(y=y, X=X))

gnames = table(stringi::stri_sub(gsub("^X.", "", colnames(X)) , 1, -2))
groups = rep(1:length(gnames), gnames)

# package all the data
d = list(y=y, X=X, groups=groups, gnames=gnames, snps=snps, genotypes=as.matrix(genotypes))

rm(list=c("i", "p", "X", "S", "high_cor", "groups","ids", "keep", "y", "varXq25", 
	  "snps", "mafXq25", "mafX","gnames", "genotypes", "phenotypes"))


# -------------------------------------------------------------------------------
# 			Create Train and Test sets
# -------------------------------------------------------------------------------
# get xval indices
xval = cv(length(d$y), 1)

d.train = list(y=d$y[xval$train], X=d$X[xval$train, ], groups=d$groups, 
	       genotypes=d$genotypes[xval$train, ])
d.train = std(d.train)

d.test  = list(y=d$y[xval$test ], X=d$X[xval$test , ], groups=d$groups, 
	       genotypes=d$genotypes[xval$train, ])


# -------------------------------------------------------------------------------
# 			    Fit the models
# -------------------------------------------------------------------------------
f1 = gsvb::gsvb.fit(d.train$y., d.train$X., d$groups, intercept=F, niter=20,
		    track_elbo=FALSE)
f2 = gsvb::gsvb.fit(d.train$y., d.train$X., d$groups, intercept=F, 
		    diag_covariance=FALSE, niter=30, track_elbo=FALSE)

save(list=c("f1", "f2"), file="../../rdata/application/mouse/models.RData")

layout(matrix(1:2, nrow=1))
plot(f1$g); plot(f2$g)
plot(f1$beta); plot(f2$beta)

f1$beta[f1$g[d$groups] > 0.1]

# which groups (SNPs) are selected
which(f1$g > 0.5)
which(f2$g > 0.5)

get_genes(d$snps[f1$g > 0.1])
get_genes(d$snps[f2$g > 0.1])

snps[which(f1$g > 0.1)]
snps[which(f2$g > 0.1)]

mse(d.train$y., d.train$X. %*% f1$beta)
mse(d.train$y., d.train$X. %*% f2$beta)


# -------------------------------------------------------------------------------
# 			    Evaluate the models
# -------------------------------------------------------------------------------
# rescale so they work on the test data
f1. = rescale_fit(f1, d.train)
f2. = rescale_fit(f2, d.train)

round(f1.$beta[f1$g[d$groups] > 0.1], 4)
round(f2.$beta[f2$g[d$groups] > 0.1], 4)

sum(f1$g > 0.1)
sum(f2$g > 0.1)

mse(d.test$y, d.test$X %*% f1.$beta + f1.$intercept)
mse(d.test$y, d.test$X %*% f2.$beta + f2.$intercept)

set.seed(1234)
y.star1 = pred.y(f1., d.test$X, d.train$y, d.train$X, sample=1e4)
pp = t(apply(y.star1, 1, quantile, probs=c(0.025, 0.975)))
mean(d.test$y >= pp[ , 1] & d.test$y <= pp[, 2])
mean(pp[ , 2] - pp[ , 1])

y.star2 = pred.y(f2., d.test$X, d.train$y, d.train$X, sample=1e4)
pp = t(apply(y.star2, 1, quantile, probs=c(0.025, 0.975)))
mean(d.test$y >= pp[ , 1] & d.test$y <= pp[, 2])
mean(pp[ , 2] - pp[ , 1])

# -------------------------------------------------------------------------------
# 			   Plot the uncertainty
# -------------------------------------------------------------------------------
pdf("../figs/mouse.pdf", width=12, height=5)

strt = 30
layout(matrix(rep(1:8, times=c(7,6,6,6,7,6,6,6)), byrow=T, nrow=2))
# layout(matrix(rep(1:8, times=c(4,3,3,3,4,3,3,3)), byrow=T, nrow=2))
# layout(matrix(1:8, byrow=T, nrow=2))
for (i in 1:8) 
{
    if (i <= 4) {
	gt4 = FALSE
	y.star = y.star1[strt + i, ]
	color = colorRampPalette(c("white", "darkgreen"))(4)
    } else {
	gt4 = TRUE
	i = i - 4
	y.star = y.star2[strt + i, ]
	color = colorRampPalette(c("white", "darkblue"))(4)
    }

    denY = density(y.star)
    qs = quantile(ecdf(y.star), probs=c(0.005, 0.025, 0.05, 0.95, 0.975, 0.995))
    
    if (i == 1) {
	par(mar=c(3, 4, 1, 1.0), family="Times")
    } else {
	par(mar=c(3, 2, 1, 1.0), family="Times")
    }

    plot(denY,  ylim=c(0.15, 4.5), xlim=c(0.0, 1.0), yaxt="n", xaxt="n",
	 main="", lwd=0, ylab="", xlab="")
    axis(2, las=2, lwd=0, lwd.ticks=0.4)
    axis(1, lwd=0, lwd.ticks=0.4)
    grid(5, 5)
    mtext("y", 1, line=1.4, cex=0.7)

    if (i == 1) {
	mtext("density", 2, line=1.5, cex=0.7)
	mtext(ifelse(gt4, "GSVB-B", "GSVB-D"), 2, line=2.7, font=2)
    }

    mean.y = mean(y.star)
    mean.y.yval = denY$y[sum(denY$x <= mean.y)]
    lines(c(mean.y, mean.y), c(0, mean.y.yval), col=adjustcolor(color[4], .7), lwd=3)

    lines(denY, lwd=2)
    color = adjustcolor(color, 0.3)

    xx = denY$x[denY$x > qs[3] & denY$x < qs[4]]
    yy = denY$y[denY$x > qs[3] & denY$x < qs[4]]
    polygon(c(xx, rev(xx)), c(yy, rep(0, length(yy))), col=color[4], border=F)

    xx = denY$x[denY$x > qs[2] & denY$x < qs[3]]
    yy = denY$y[denY$x > qs[2] & denY$x < qs[3]]
    polygon(c(xx, rev(xx)), c(yy, rep(0, length(yy))), col=color[3], border=F)

    xx = denY$x[denY$x > qs[4] & denY$x < qs[5]]
    yy = denY$y[denY$x > qs[4] & denY$x < qs[5]]
    polygon(c(xx, rev(xx)), c(yy, rep(0, length(yy))), col=color[3], border=F)

    xx = denY$x[denY$x > qs[1] & denY$x < qs[2]]
    yy = denY$y[denY$x > qs[1] & denY$x < qs[2]]
    polygon(c(xx, rev(xx)), c(yy, rep(0, length(yy))), col=color[2], border=F)

    xx = denY$x[denY$x > qs[5] & denY$x < qs[6]]
    yy = denY$y[denY$x > qs[5] & denY$x < qs[6]]
    polygon(c(xx, rev(xx)), c(yy, rep(0, length(yy))), col=color[2], border=F)

    abline(v = d.test$y[strt + i], col=2, lwd=3)
    text(0.90, y=4.3, labels=paste0("y*: ", d.test$y[strt + i]))
    text(0.83, y=4.0, labels=expression(hat(y)))
    text(0.93, y=4.0, label=sprintf(": %.2f", mean.y))
    
    if (!gt4) {
	mtext(paste0(LETTERS[i], ".1"), 2, las=1, padj=-8.2, adj=1.3, font=2)
    } else {
	mtext(paste0(LETTERS[i], ".2"), 2, las=1, padj=-8.2, adj=1.3, font=2)
    }
}
dev.off()

as.hexmode(col2rgb(adjustcolor(colorRampPalette(c("white", "darkblue"))(4), 0.3)))


# -------------------------------------------------------------------------------
# 			   	Misc.
# -------------------------------------------------------------------------------
lambda0 = seq(50, 5, by=-5)
cv.ssgl = sparseGAM::cv.SSGL(d.train$y, d.train$X, d.train$groups, lambda0=lambda0)

f.ssgl = sparseGAM::SSGL(d.train$y, d.train$X, d.train$X, 
			 d.train$groups, lambda0=5)

plot(f.ssgl$beta)
sum(nzero.groups(f.ssgl$beta, d$groups))
mse(d.test$y, d.test$X %*% f.ssgl$beta + f.ssgl$beta0)


save(list=c("cv.ssgl", "f.ssgl"), file="../../rdata/application/mouse/models_ssgl.RData")

# -------------------------------------------------------------------------------
# 			   	Misc.
# -------------------------------------------------------------------------------
# fit as univariate analysis
# for this we are going to use the spike and slab prior

# fit the method
ff = sparsevb::svb.fit(d.train$genotypes, d.train$y., family="linear")

# compute the posterior mean and the intercept
ff$beta = ff$mu * ff$g
ff$intercept = mean(d.train$y - d.train$genotypes %*% ff$beta)

# compute the held out MSE
mse(d.test$y, d.test$genotypes %*% ff$beta + ff$intercept)
mse(d.test$y, d.test$genotypes %*% ff$beta + ff$intercept)
get_genes(d$snps[ff$g > 0.1])


save(list=c("ff"), file="../../rdata/application/mouse/models_svb.RData")
