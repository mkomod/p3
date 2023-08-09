rm(list=ls())

pheno = fread("../../data/breast_cancer/TCGA-BRCA.GDC_phenotype.tsv")
genex = fread("../../data/breast_cancer/TCGA-BRCA.htseq_fpkm.tsv")


# get the gene names
genes = genex[, 1]
genex = genex[ , -1]

# get the patient IDs
genex_ids = colnames(genex)
length(genex_ids)

# get the patient ids for the phenotype data
pheno_ids = pheno$submitter_id.samples
pheno_ids = pheno_ids[pheno$new_tumor_event != ""]
pheno = pheno[pheno$new_tumor_event != "", ]

# extract the ids that are shared in both datasets
id_intersection = intersect(pheno_ids, genex_ids)

# subset the datasets
ss = pheno_ids %in% id_intersection
y = pheno[ss, ]
y = y$new_tumor_event
y = y[order(pheno_ids[ss])]
y = as.numeric(y == "YES")

ss = which(genex_ids %in% id_intersection)
X = genex[ , ..ss]
X = t(X)
X = X[order(genex_ids[ss]), ]


# keep 5000 genes with highest variance (on the log scale)
X_var = apply(X, 2, var)
keep = order(X_var, decreasing = T)[1:1000]
X = X[ , keep]


length(y)
dim(X)

G = ncol(X)
n =  nrow(X)
mg = 3

Xnew = matrix(NA, nrow=n, ncol=G*mg)

for (g in 1 : G) {
    splineTemp = ns(X[,g], df=mg)
    for (m in 1 : mg) {
	Xnew[,mg*(g-1) + m] = splineTemp[,m]
    }
}

groups = rep(1:G, each=mg)

X = scale(X, center=T)
f = gsvb::gsvb.fit(y, X, 1:1000, family="binomial-jaakkola", lambda=0.03, intercept=FALSE, track_elbo=FALSE)
plot(f$g)


plot(Xnew[ , groups == 335] %*% f$beta_hat[groups == 335], y)


1:10 %*% ns(1:10, df=3) 

ns(1:100, df=3)

x = seq(-1, 1, length.out=20)
y = x^3 + rnorm(length(x), 0, 0.1)
plot(x, y)

f = lm(y ~ ns(x, df=5))
dn =data.frame(x = seq(-1, 1, length.out=50))

plot(x, y)
lines(unlist(dn), predict(f, newdata=dn))
