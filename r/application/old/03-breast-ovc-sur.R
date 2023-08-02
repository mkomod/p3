library(data.table)
library(splines)


br_gx


ov_gx = data.table::fread("../../data/ovarian_cancer/tcga_ov_affy.tsv")
ov_pheno = data.table::fread("../../data/ovarian_cancer/tcga_ov_pheno.tsv")




ov_pheno_id = ov_pheno$sampleID
ov_pheno_li = ov_pheno$lymphatic

to_keep = (ov_pheno_li != "")

ov_pheno_id = ov_pheno_id[to_keep]
ov_pheno_li = ov_pheno_li[to_keep]

ov_gx_id = colnames(ov_gx)

shared_id = intersect(ov_gx_id, ov_pheno_id)

to_keep = ov_gx_id %in% shared_id
ov_gx_id = ov_gx_id[to_keep]
ov_gx = ov_gx[ , ..to_keep]

data.table::setcolorder(ov_gx, sort(ov_gx_id))

to_keep = ov_pheno_id %in% shared_id
ov_pheno_id = ov_pheno_id[to_keep]
ov_pheno_li = ov_pheno_li[to_keep]
ov_pheno_li = ov_pheno_li[order(ov_pheno_id)]

sort(ov_gx_id) == sort(ov_pheno_id)


X1 = data.table::copy(ov_gx)
y1 = ov_pheno_li

y1 = as.numeric(y1 == "YES")
table(y1)

X1 = as.matrix(X1)
X1 = t(X1)

to_rm = sample(which(y1 == 1), 57)
y1 = y1[-to_rm]
X1 = X1[-to_rm, ]

varX1 = apply(X1, 2, var)
to_keep = varX1 > 0.5
X1 = X1[ , to_keep]
X1 = scale(X1)

Z = matrix(nrow=nrow(X1), ncol=0)
d = 4
for (j in 1:ncol(X1)) {
    Z = cbind(Z, splines::bs(X1[ ,j], df=d))
}
Z = scale(Z, center=T, scale=F)
groups = rep(1:ncol(X1), each=d)


ff = gsvb::gsvb.fit(y1, Z, groups, family="binomial-jaak", diag_covariance=TRUE,
		    niter=5, track_elbo=FALSE)

ff
plot(ff$g)
which.max(ff$g)

plot(ff$beta_hat)
plot(ff$mu)

cor(Z[ , 1:4])
