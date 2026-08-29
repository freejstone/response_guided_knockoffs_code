###########
#simulation functions
###########
#erdos_renyi_matrix: Function for constructing design matrix with ER covariance matrix. (based off knockpy.dgp.ErdosRenyi)
#ER_cov: Converts output of erdos_renyi_matrix to a correlation matrix. (based off lines from knockpy.dgp.DGP.sample_data) 
#AR1: Function for constructing design matrix with AR1 covariance matrix (based off knockpy.dgp.AR1)
#create_sparse_beta_contiguous: Creates a true beta coefficients for linear model. (based off knockpy.dgp.create_sparse_coefficients)
erdos_renyi_matrix = function(p,
                              delta = 0.2,         
                              lower = 0.1, upper = 1,
                              tol = 1e-1) {
  V = matrix(0, p, p)
  idx = upper.tri(V, diag = TRUE)
  triang_size = sum(idx)
  
  mask  = rbinom(triang_size, size = 1, prob = delta)
  vals  = runif(triang_size, min = lower, max = upper)
  flags = sample(c(-1, 1), size = triang_size, replace = TRUE)
  
  triang = mask * flags * vals
  V[idx] = triang
  V = V + t(V)                         # symmetrize
  diag(V) = 1                          # set diagonal to 1
  V = shift_until_PSD(V, tol = tol)    # ensure PSD (may move diag off 1)
  V
}
ER_cov = function(p,
                  delta = 0.2,        
                  lower = 0.1, upper = 1,
                  tol = 1e-1,
                  max_corr = 1) {
  V = erdos_renyi_matrix(p, delta = delta, lower = lower, upper = upper, tol = tol)
  # Convert to correlation (like knockpy.utilities.cov2corr)
  Sigma = stats::cov2cor(V)
  # Clamp off-diagonals to [-max_corr, +max_corr]
  off = row(Sigma) != col(Sigma)
  Sigma[off] = pmax(pmin(Sigma[off], max_corr), -max_corr)
  diag(Sigma) = 1
  Sigma
}
ER_prec = function(p, delta = 0.2, lower = 0.1, upper = 1,
                                 tol = 1e-1) {
  Q  = erdos_renyi_matrix(p, delta = delta, lower = lower, upper = upper, tol = tol)
  V  = solve(Q)
  stats::cov2cor(V)
}
AR1 = function(p = 30,
               a = 1,
               b = 1,
               tol = NULL,
               max_corr = 1,
               rho = NULL) {
  # Generate log rhos
  if (is.null(rho)) {
    # Draw successive correlations from Beta(a, b), truncated at max_corr
    rhos = pmin(stats::rbeta(p, shape1 = a, shape2 = b), max_corr)
    rhos_log = log(rhos)
  } else {
    if (abs(rho) > 1) {
      stop(sprintf("rho %s must be a correlation between -1 and 1", rho))
    }
    rhos_log = rep(log(rho), p)
  }
  # First entry is set to 0 (i.e. log 1)
  rhos_log[1] = 0
  
  # Cumulative logs
  cumrhos = cumsum(rhos_log)
  
  # Log correlations: -|cumrhos_i - cumrhos_j|
  log_corrs = -abs(outer(cumrhos, cumrhos, "-"))
  
  # Correlation matrix
  corr_matrix = exp(log_corrs)
  
  # Ensure PSD-ness if tol is provided
  if (!is.null(tol)) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required when 'tol' is not NULL.")
    }
    # nearPD tries to find the nearest positive definite correlation matrix
    # nearPD with corr=TRUE already returns a correlation matrix
    npd = Matrix::nearPD(corr_matrix,
                         corr    = TRUE,
                         keepDiag = TRUE,
                         eig.tol = tol)
    corr_matrix = as.matrix(npd$mat)
  }
  corr_matrix
}
block_equi_corr = function(p, block_size = 5, rho = 0.5, gamma = 0) {
  stopifnot(p %% block_size == 0)
  num_groups = p / block_size
  groups0 = rep(0:(num_groups - 1), each = block_size)  # 0-based like Python
  get_corr = function(g1, g2) ifelse(g1 == g2, rho, gamma * rho)
  G = outer(groups0, groups0, Vectorize(get_corr))
  G = G + diag(1, p) - diag(diag(G), p)
  G
}
create_sparse_coefficients = function(
    p,
    sparsity = 0.5,
    groups = NULL,
    coeff_size = 1,
    coeff_dist = c("none", "normal", "uniform", "dsliu2020", "gmliu2019"),
    sign_prob = 0.5,          # PROBABILITY OF NEGATIVE sign (matches code)
    iid_signs = TRUE,
    corr_signals = FALSE,     # if TRUE and groups=NULL: contiguous block
    n = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  coeff_dist = match.arg(coeff_dist)
  
  # Choose non-zeros
  if (!is.null(groups)) {
    ug = sort(unique(groups))
    num_groups = length(ug)
    num_nonzero_groups = floor(sparsity * num_groups)
    chosen = sample(ug, num_nonzero_groups, replace = FALSE)
    beta = ifelse(groups %in% chosen, coeff_size, 0)
  } else {
    num_nonzero = floor(sparsity * p)
    if (corr_signals) {
      start = sample.int(p - num_nonzero + 1, 1) - 1 # 0-based idea
      beta = rep(0, p)
      beta[(start + 1):(start + num_nonzero)] = coeff_size
    } else {
      beta = c(rep(coeff_size, num_nonzero), rep(0, p - num_nonzero))
      beta = sample(beta)  # shuffle
    }
  }
  
  beta_nz = beta != 0
  k = sum(beta_nz)
  
  # Signs: signs = 1 - 2 * Bernoulli(sign_prob)  => negative with prob sign_prob
  if (iid_signs) {
    signs = 1 - 2 * rbinom(p, 1, sign_prob)
    beta = beta * signs
  } else {
    num_pos = floor(k * (1 - sign_prob))
    signs = c(rep(1, num_pos), rep(-1, k - num_pos))
    signs = sample(signs)
    beta[beta_nz] = abs(beta[beta_nz]) * signs
  }
  
  # Magnitude distributions
  if (coeff_dist == "normal") {
    # sqrt(coeff_size) * N(0,1) on non-zeros
    tmp = rep(0, p)
    tmp[beta_nz] = sqrt(coeff_size) * rnorm(k)
    beta = tmp
  } else if (coeff_dist == "uniform") {
    # beta = beta * runif(p)/2 + beta/2  => |beta| ~ Unif(δ/2, δ), with sign from above
    beta = beta * runif(p)/2 + beta/2
  } else if (coeff_dist == "dsliu2020") {
    if (k != 50) stop("For dsliu2020, need exactly 50 non-zeros (per knockpy).")
    if (is.null(n)) stop("dsliu2020 requires n.")
    variance = 10 * sqrt(log(p) / n)
    tmp = rep(0, p)
    tmp[beta_nz] = sqrt(variance) * rnorm(k)
    beta = tmp
  } else if (coeff_dist == "gmliu2019") {
    if (k != 60) stop("For gmliu2019, need exactly 60 non-zeros (per knockpy).")
    if (is.null(n)) stop("gmliu2019 requires n.")
    variance = 20 / sqrt(n)
    tmp = rep(0, p)
    tmp[beta_nz] = sqrt(variance) * rnorm(k)
    beta = tmp
  } # "none": leave magnitudes at coeff_size (with sign)
  
  beta
}
