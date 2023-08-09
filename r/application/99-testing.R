# ------------------------------------------------------------------------------------------
# 					Testing
# ------------------------------------------------------------------------------------------
genx = fread("../../data/breast_cancer/TCGA-BRCA.htseq_fpkm.tsv")

genx

R = cor(X)

proc_R = process_LD(R, thresh=0.4)


process_LD <- function(LD, thresh=0.95)
{
    variant_names = colnames(LD)
    
    p = ncol(LD)
    
    selected = rep(FALSE, p)
    ld_pairs = list()
    ld_pairs[p + 1] = NULL
    
    for (i in 1:(nrow(LD) - 1))
    {
        if (!selected[i])
        {
            # which elements are in LD with i
            in_ld_with_i = which(LD[i, (i+1):p] >= thresh) + i
            ld_pairs[[i]] = in_ld_with_i
            
            # elements which are in LD with i are then considered to be the
            # same as i
            selected[in_ld_with_i] <- TRUE
        }
    }
    
    if (!is.null(variant_names)) {
        names(selected) = variant_names
        names(ld_pairs) = variant_names
    }
    
    # TODO: these should be re-named
    return(list(
        main_variants = which(!selected),
        in_ld_with_main_variants = ld_pairs
    ))
}


n = 100
p = 10
X = matrix(rnorm(n * p, 0, 1), nrow=n)
y = cos(1.5 * X[ , 2]) + rnorm(n, 0, 0.2)
plot(y)




df = 4
Xnew = matrix(nrow=n, ncol=p*df)
for (j in 1:p) {
    Xnew[ , (df*(j-1)+1):(j*df)] = ns(X[ , j], df=df)
}
groups = rep(1:10, each=df)

f = gsvb::gsvb.fit(y, Xnew, groups, intercept=T, niter=350, 
		   diag_covariance=FALSE)
plot(f$g)
plot(f$beta_hat)

xs = seq(min(X[ , 2]), max(X[ , 2]), length.out=500)
newd = cbind(matrix(0, nrow=500, ncol=df),
	     ns(xs, df=df),
	     matrix(0, nrow=500, ncol=df * (p-2)))
pred = gsvb::gsvb.predict(f, newd)

plot(xs, pred$mean, type="l", lwd=3, ylim=c(-3, 4))
lines(xs, pred$quantile[1, ])
lines(xs, pred$quantile[2, ])
points(X[ , 2], y)


