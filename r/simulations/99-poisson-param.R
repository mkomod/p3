library(spsl)
library(ParBayesianOptimization)
library(coda)
library(parallel)
source("./00-functions.R")


objective_function <- function(k1, k2) {
  N = 10
  scores <- numeric(N)
  for (i in 1:N) {
    tryCatch({
        seed = sample(1:100, 1)
        d <- dgp_block(400, 1000, 5, 2, 0.45, list(model="poisson", corr=0.6, block_size=50), seed=seed)
        m_par = list(family="poisson", lambda=1, a0=1, b0=1000/5 + 1, a_t=1e-3, b_t=1e-3,
               mcmc_samples=2e4, burnin=5e3, intercept=FALSE, kp_1=k1, kp_2=k2)
        f <- m_spsl(d, m_par)
        scores[i] <- -f[length(f) - 1]
    }, error = function(e) {
        cat("Error in iteration", i, ":", e$message, "\n")
        scores[i] <- NA
    })
  }
  return(mean(scores, na.rm=TRUE))
}

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
k1_values <- seq(0.008, 0.020, by = 0.002)
k2_values <- seq(24, 32, by = 1)
length(k1_values) * length(k2_values)

# Perform the grid search
grid_search_results <- grid_search(k1_values, k2_values)

print(grid_search_results)
sapply(grid_search_results$AllParams, function(x) x$score) > -2

which.max(sapply(grid_search_results$AllParams, function(x) x$score) )
expand.grid(k1_values, k2_values)[22, ]
grid_search_results$AllParams[[57]]

# Create a matrix to store the scores
score_matrix <- matrix(NA, nrow = length(k1_values), ncol = length(k2_values))
rownames(score_matrix) <- k1_values
colnames(score_matrix) <- k2_values

# Fill the matrix with the scores from the grid search results
for (result in grid_search_results$AllParams) {
  k1_idx <- which(k1_values == result$k1)
  k2_idx <- which(k2_values == result$k2)
  score_matrix[k1_idx, k2_idx] <- result$score
}

# Print the score matrix
print(score_matrix)
rowMeans(score_matrix)
colMeans(score_matrix)

0.022 / 22
0.02 / 14
0.024 // 12
0.02 // 16


