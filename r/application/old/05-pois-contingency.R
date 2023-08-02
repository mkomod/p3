site <- gl(n = 3, k = 1, length = 12) 
# gl generates levels of a factor
tumor <- gl(n = 4, k = 3) #each 3
cases <- c(22, 2, 10, 16, 54, 115, 19, 33, 73, 11, 17, 28)
cancer <- data.frame(tumor, site, cases)

cancer

cancer.m0 <- glm(cases ~ 1, family = poisson, data = cancer)
cancer.m1 <- glm(cases ~ tumor, family = poisson, data = cancer)
cancer.m2 <- glm(cases ~ site, family = poisson, data = cancer)
cancer.m3 <- glm(cases ~ tumor + site, family = poisson, data = cancer)
# Saturated model
cancer.m4 <- glm(cases ~ tumor * site, family = poisson, data = cancer)
# Analysis of deviance
# Same syntax as for GLM
drop1(cancer.m4)

summary(cancer.m4)

cancer.m4$model


tumor == 1
site == 1



sum(cases[site == 1])
sum(cases[site == 2])
sum(cases[site == 3])

X = model.matrix(cancer.m4)


gsvb::gsvb.fit(cases, 
