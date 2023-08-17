# -------------------------------------------------------------------------------
# 			     Create the groups
# -------------------------------------------------------------------------------
# keep the snps which are in the snp database
keep = snps %in% snp_locations$marker
snps = snps[keep]
snp_locations = snp_locations[ snp_locations$marker %in% snps, ]

# the order is the same
all(snp_locations$marker == snps) 

# find the Genes associated with each SNP
marker_locations = marker_locations[marker_locations$Marker.Type == "Gene", ]

genes = apply(snp_locations, 1, function(snp) 
{
    # make sure the chromosome is the same
    marker_locations_sub = marker_locations[marker_locations$Chrom == snp[3], ]

    # is the location of the SNP wiithin the gene location
    bp = as.numeric(snp[4])
    marker_locations_sub$Marker.Symbol[
	which(marker_locations_sub$Start < bp & bp < marker_locations_sub$End)
    ]
})

# there are multiple genes to each snp
# first we remove the snps that have no genes
keep = unlist(lapply(genes, function(g) length(g) > 0))

genes = genes[keep]
snps_locations = snp_locations[keep, ]
snps = snps[keep]

# some snps are associated to more than one gene 
# we remove the least common gene
tg = table(unlist(genes))
# tg = tg[which(tg > 1)]
genes = lapply(genes, function(g) {
    g[g == names(which.max(tg[names(tg) %in% g]))]
})

# remove genes with no snps
keep = unlist(lapply(genes, function(g) length(g) > 0))
genes = genes[keep]
snps_locations = snp_locations[keep, ]
snps = snps[keep]


# order the design by the genes
o = order(unlist(genes))
genes = genes[o]
genes = unlist(genes)
snps = snps[o]


# remove snps in the design with no group
genotypes = genotypes[ , colnames(genotypes) %in% snps]
genotypes = genotypes[ , snps]


# create the groups
X = as.matrix(genotypes)
X = apply(X, 2, as.factor)
XX = model.matrix(y ~ . - 1, data.frame(y=y, X=X))

groups = as.numeric(factor(genes))
groups = rep(groups, table(stringi::stri_sub(gsub("^X.", "", colnames(XX)), 1, -2)))

# package the data
d = list(y = y - mean(y), X=XX, groups = groups, intercept=mean(y))

d = std(d)

# -------------------------------------------------------------------------------
# 			     	Fit the methods
# -------------------------------------------------------------------------------
f1 = gsvb::gsvb.fit(d$y, d$X, d$groups, intercept=F, niter=20, track_elbo=FALSE)
f2 = gsvb::gsvb.fit(d$y, d$X, d$groups, diag_covariance=FALSE, intercept=T, 
		    niter=20, track_elbo=FALSE)
plot(f1$g)
plot(f1$beta)

plot(f2$g)
plot(f2$beta)

unique(genes)[which(f1$g > 0.5)]













# -------------------------------------------------------------------------------
# 			  Use genes to create groups
#
# LATEST 07.08.2023
# -------------------------------------------------------------------------------
#
#
#
# keep the snps which are in the snp database
keep = snps %in% snp_locations$marker
snps = snps[keep]
snp_locations = snp_locations[ snp_locations$marker %in% snps, ]

# the order is the same
all(snp_locations$marker == snps) 

# find the Genes associated with each SNP
marker_locations = marker_locations[marker_locations$Marker.Type == "Gene", ]

genes = apply(snp_locations, 1, function(snp) 
{
    # make sure the chromosome is the same
    marker_locations_sub = marker_locations[marker_locations$Chrom == snp[3], ]

    # is the location of the SNP wiithin the gene location
    bp = as.numeric(snp[4])
    marker_locations_sub$Marker.Symbol[
	which(marker_locations_sub$Start < bp & bp < marker_locations_sub$End)
    ]
})

# there are multiple genes to each snp
# first we remove the snps that have no genes
keep = unlist(lapply(genes, function(g) length(g) > 0))

genes = genes[keep]
snps_locations = snp_locations[keep, ]
snps = snps[keep]

# some snps are associated to more than one gene 
# we remove the least common gene
tg = table(unlist(genes))
# tg = tg[which(tg > 1)]
genes = lapply(genes, function(g) {
    g[g == names(which.max(tg[names(tg) %in% g]))]
})

# remove genes with no snps
keep = unlist(lapply(genes, function(g) length(g) > 0))
genes = genes[keep]
snps_locations = snp_locations[keep, ]
snps = snps[keep]


# order the design by the genes
o = order(unlist(genes))
genes = genes[o]
genes = unlist(genes)
snps = snps[o]


# remove snps in the design with no group
genotypes = genotypes[ , colnames(genotypes) %in% snps]
genotypes = genotypes[ , snps]


# create the groups
X = as.matrix(genotypes)
groups = as.numeric(factor(genes))

# package the data
d = list(y=y, X=X, groups=groups)
d = std(d)

# -------------------------------------------------------------------------------
# 			     	Fit the methods
# -------------------------------------------------------------------------------
f1 = gsvb::gsvb.fit(d$y., d$X., d$groups, intercept=F, niter=100, track_elbo=FALSE)
f2 = gsvb::gsvb.fit(d$y, d$X, d$groups, diag_covariance=FALSE, intercept=T,
		    niter=50, track_elbo=FALSE)

layout(matrix(1:2, nrow=1))
plot(f1$g)
plot(f2$g)

f1 = gsvb::gsvb.fit(d$y, d$X, d$groups, intercept=F, niter=20, track_elbo=FALSE)
f2 = gsvb::gsvb.fit(d$y, d$X, d$groups, diag_covariance=FALSE, intercept=T, 
		    niter=20, track_elbo=FALSE)
plot(f1$g)
plot(f1$beta)

plot(f2$g)
plot(f2$beta)

unique(genes)[which(f1$g > 0.5)]
unique(genes)[which(f2$g > 0.5) - 1]




