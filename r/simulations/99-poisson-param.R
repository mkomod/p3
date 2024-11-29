library(spsl)
library(ParBayesianOptimization)
library(coda)

# setwd("./r/simulations")
# source("00-functions.R")
# 
# objective_function <- function(k1, k2) {
#   seed = sample(1:100, 1)
# 
#   d <- dgp_diag(400, 1000, 5, 2, 0.45, list(model="poisson", corr=0), seed=seed)
#   m_par = list(family="poisson", lambda=1, a0=1, b0=1000/5 + 1, a_t=1e-3, b_t=1e-3,
#          mcmc_samples=1e4, burnin=1e3, intercept=TRUE, kp_1=0.10, kp_2=25.0)
#   f <- m_spsl(d, m_par)
#   f[length(f) - 1]
# 
# 
#   
#   # Fit the model with the given kernel parameters
#   f0 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
#   f1 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
#   f2 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
#   f3 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
#   
#   # Create MCMC list
# 	l =	coda::mcmc.list(
#     coda::mcmc(t(
# 			f0$B * f0$Z[f0$parameters$groups, ]	
# 		)),
# 		coda::mcmc(t(
# 			f1$B * f1$Z[f1$parameters$groups, ]	
# 		)),
# 		coda::mcmc(t(
# 			f2$B * f2$Z[f2$parameters$groups, ]	
# 		)),
# 		coda::mcmc(t(
# 			f3$B * f3$Z[f3$parameters$groups, ]	
#     ))
#   )
# 
#   
#   # Calculate Gelman-Rubin diagnostic
#   tryCatch({
#     gg <- coda::gelman.diag(l, multivariate=FALSE)
#     max_psrf = max(gg$psrf[ , 1], na.rm=TRUE)
# 
#     # Return the multivariate potential scale reduction factor (mpsrf) as the objective
#     # return(list(Score = -gg$mpsrf))  # Negative because we want to minimize mpsrf
#     return(list(Score = -max_psrf))
#   }, error = function(e) {
#     cat("Error in calculating Gelman-Rubin diagnostic: ", e$message, "\n")
#     return(list(Score = -500))
#   })
# }

library(spsl)
library(ParBayesianOptimization)
library(coda)



objective_function <- function(k1, k2) {
  seed = sample(1:100, 1)
  d <- dgp_diag(400, 1000, 5, 2, 0.45, list(model="poisson", corr=0), seed=seed)
  m_par = list(family="poisson", lambda=1, a0=1, b0=1000/5 + 1, a_t=1e-3, b_t=1e-3,
         mcmc_samples=1e4, burnin=1e3, intercept=TRUE, kp_1=k1, kp_2=k2)
  f <- m_spsl(d, m_par)
  return(list(Score = -f[length(f) - 1]))
}


opt_results_pois <- bayesOpt(
  FUN = objective_function,
  bounds = list(k1 = c(0.05, 0.15), k2 = c(10, 30)),
  initPoints = 5,
  iters.n = 5,
  acq = "ei"
)

# Print the optimal parameters
opt_results_pois$scoreSummary
ParBayesianOptimization::getBestPars(opt_results_pois)
# objective_function(0.13, 23)
