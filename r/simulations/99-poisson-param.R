library(spsl)
library(ParBayesianOptimization)
library(coda)

source("00-functions.R")

objective_function <- function(k1, k2) {
  seed = sample(1:100, 1)
  d <- dgp_diag(400, 1000, 5, 2, bmax=0.45, list(model="poisson", corr=0.6), seed=seed)
  
  # Fit the model with the given kernel parameters
  f0 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  f1 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  f2 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  f3 <- spsl::spsl.fit(d$y, d$X, d$groups, family="poisson", mcmc_samples=1.5e4, burnin=1e3, kernel_param_1=k1, kernel_param_2=k2)
  
  # Create MCMC list
  l <- coda::mcmc.list(
    coda::mcmc(t(f0$B)),
    coda::mcmc(t(f1$B)),
    coda::mcmc(t(f2$B)),
    coda::mcmc(t(f3$B))
  )
  
  # Calculate Gelman-Rubin diagnostic
  tryCatch({
    gg <- coda::gelman.diag(l, multivariate=TRUE)
    
    # Return the multivariate potential scale reduction factor (mpsrf) as the objective
    return(list(Score = -gg$mpsrf))  # Negative because we want to minimize mpsrf
  }, error = function(e) {
    cat("Error in calculating Gelman-Rubin diagnostic: ", e$message, "\n")
    return(list(Score = -500))
  })
}

opt_results_pois <- bayesOpt(
  FUN = objective_function,
  bounds = list(k1 = c(0.05, 0.15), k2 = c(17, 30)),
  initPoints = 5,
  iters.n = 20,
  acq = "ei"
)

# Print the optimal parameters
opt_results_pois$scoreSummary
ParBayesianOptimization::getBestPars(opt_results_pois)

