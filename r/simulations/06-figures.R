

metric <- c("coverage.non_zero", "length.non_zero", "coverage")
metric.title <- c("Coverage Non-zero", "Length Non-zero", "PP Coverage")
for (i in 1:3) 
{
    # pdf(sprintf("../figs/cov_comp_gaus_%d.pdf", i), width=12, height=3)
    # png(sprintf("../figs/cov_comp_gaus_%d.png", i), width=12, height=3, units="in", res=200)
    layout(matrix(1:3, nrow=1, byrow=T))
    gsize <- g[i]; n0g <- s[i]; nn <- n[i]
    for (m in 1:3) 
    {
	met <- metric[m]
	ddat <- dat[ dat[ , "g"] == gsize & dat[ , "s"] == n0g 
		    & dat[ , "n"] == nn, ]
	f <- as.formula(sprintf("%s ~ m*d", met))

	par(mar=c(3, 3, 3, 0))
	boxplot(f, data=ddat, main=metric.title[m], boxwex=0.2, col=c(2,3))
    }
    # dev.off()
}



metrics <- c("l2", "auc")
metric.title <- c("l2-error", "AUC")
simnums <- 1:2
layout(matrix(1:(length(metrics) * length(simnums)), 
	      ncol=length(metrics), byrow=T))
for (sim in simnums)
{
    nn <- n[sim]; pp <- p[sim]; gg <- g[sim]; ss <- s[sim]
    for (met in metrics)
    {
	mm <- metrics[met]
	ddat <- dat[  dat[ , "g"] == gg & dat[ , "s"] == ss 
		    & dat[ , "n"] == nn & dat[ , "p"] == pp, ]

	f <- as.formula(sprintf("%s ~ m*d", met))

	par(mar=c(3, 3, 3, 0))
	boxplot(f, data=ddat, boxwex=0.2, col=1:3)
	grid()
    }
}

library(ggplot2)  # Load ggplot2 package

# Create the plot using ggplot2
pplot <- ggplot2::ggplot(data.frame(dat), aes(x=methods)) +
  geom_histogram(binwidth = 2, color = "black", fill = "skyblue", alpha = 0.8) +
  labs(title = "Histogram of Data", x = "Values", y = "Frequency") +
  theme_minimal()

# Display the plot
print(plot)

