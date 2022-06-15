# Misc code for debugging
source("./00-functions.R")

d <- dgp_diag(200, 1000, 5, 3, list(corr=0))
m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
f <- m_gsvb(d)


d <- dgp_block(200, 1000, 5, 3, list(corr=0.6, block_size=50), seed=91)
m_par <- list(lambda=1, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
m_par <- list(lambda=0.5, a0=1, b0=200, a_t=1e-3, b_t=1e-3)
f <- m_gsvb(d, m_par=m_par)
