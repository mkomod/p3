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




d <- dgp_diag(400, 100, 5, 2, 0.45, list(model="poisson", corr=0.6), seed=1)
m_par = list(family="poisson", lambda=1, a0=1, b0=1000/5 + 1, a_t=1e-3, b_t=1e-3,
       mcmc_samples=2e4, burnin=1e4, intercept=FALSE, kp_1=k1, kp_2=k2)
f <- m_spsl(d, m_par)
f0 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=0.06, kernel_param_2=10)
f1 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
f2 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
f3 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
#

f0$B[ (f0$g > 0.5)[f0$parameters$groups], ] 


library(spsl)
library(ParBayesianOptimization)
library(coda)
library(parallel)

source("./00-functions.R")

objective_function <- function(k1, k2) {
  scores <- numeric(5)
  for (i in 1:5) {
    seed = sample(1:100, 1)
    d <- dgp_diag(400, 1000, 5, 2, 0.45, list(model="poisson", corr=0.6), seed=seed)
    m_par = list(family="poisson", lambda=1, a0=1, b0=1000/5 + 1, a_t=1e-3, b_t=1e-3,
           mcmc_samples=2e4, burnin=5e3, intercept=FALSE, kp_1=k1, kp_2=k2)
    f <- m_spsl(d, m_par)
    scores[i] <- -f[length(f) - 1]
  }
  return(mean(scores))
}

# objective_function(0.01, 20)


grid_search <- function(k1_values, k2_values) {
  param_grid <- expand.grid(k1 = k1_values, k2 = k2_values)
  
  scores <- mclapply(1:nrow(param_grid), function(idx) {
    k1 <- param_grid[idx, "k1"]
    k2 <- param_grid[idx, "k2"]
    sprintf("%f %f", k1, k2)
    score <- objective_function(k1, k2)
    return(list(k1 = k1, k2 = k2, score = score))
  }, mc.cores = 64)
  
  best_score <- -Inf
  best_params <- list(k1 = NA, k2 = NA)
  all_params <- list()
  
  for (result in scores) {
    all_params <- append(all_params, list(result))
    if (result$score > best_score) {
      best_score <- result$score
      best_params$k1 <- result$k1
      best_params$k2 <- result$k2
    }
  }
  
  return(list(BestScore = best_score, BestParams = best_params, AllParams = all_params))
}

# Define the range of k1 and k2 values for the grid search
k1_values <- seq(0.005, 0.05, by = 0.005)
k2_values <- seq(2, 10, by = 2)

length(k1_values) * length(k2_values)

# Perform the grid search
grid_search_results <- grid_search(k1_values, k2_values)

print(grid_search_results)
# Find the largest score that is not Inf
sapply(grid_search_results$AllParams, function(x) x$score) > -2

grid_search_results$AllParams[[1]]
