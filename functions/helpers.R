###########
#helper functions
###########
#seed_from_params: Returns a seed based off two input parameters
#save_sim_h5: Saves output file for Python
#calc_mineig: Returns smallest eigenvalue. (based off knockpy.utilities.calc_mineig)
#shift_until_PSD: Ensure diag shifts so that min eigenvalue ≥ tol. (based off knockpy.utilities.shift_until_PSD)
#scale_until_PSD: Ensure 2*Sigma - gamma S is PSD by adjusting gamma. (based off knockpy.utilities.scale_until_PSD)
#mvr_loss: Minimum-variance reconstructability objective. (based off knockpy.mrc.mvr_loss)
#maxent_loss: Objective function for entropy knockoffs. (based off knockpy.mrc.maxent_loss)

seed_from_params <- function(a, b, base_seed = 0L) {
  if (length(a) != 1L || length(b) != 1L) {
    stop("a and b must be length-1 values.")
  }
  
  # Stable string key for the pair (a, b)
  fmt <- function(x) {
    if (is.numeric(x)) sprintf("%.17g", x) else as.character(x)
  }
  key <- paste0("a=", fmt(a), "|b=", fmt(b))
  
  # Hash into 1..(2^31-1) safely (no integer overflow)
  mod <- 2147483647  # 2^31 - 1
  h <- 0             # numeric (double)
  
  raw <- charToRaw(key)
  for (x in raw) {
    h <- (h * 31 + as.integer(x)) %% mod
  }
  
  seed <- (h + (as.numeric(base_seed) %% mod)) %% mod
  seed <- as.integer(seed)
  if (is.na(seed) || seed == 0L) seed <- 1L
  seed
}
save_sim_h5 <- function(path, y, X, beta, meta = list()) {
  f <- H5File$new(path, mode = "w")
  f[["y"]]    <- as.numeric(y)
  f[["X"]]    <- as.matrix(X)
  f[["beta"]] <- as.numeric(beta)
  
  # Store metadata as HDF5 attributes (safe way)
  if (length(meta)) {
    for (nm in names(meta)) {
      h5attr(f, nm) <- meta[[nm]]
    }
  }
  f$close_all()
}
calc_mineig = function(M) {
  min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
}
shift_until_PSD = function(S, tol = 1e-5) {
  eig = eigen(S, symmetric = TRUE)
  lam = eig$values
  if (min(lam) >= tol) return(S)
  shift = tol - min(lam)
  S + diag(shift, nrow(S))
}
scale_until_PSD = function(Sigma, S, nsko = 1, tol = 1e-8, num_iter = 20) {
  lower = 0          # largest feasible gamma so far
  # This is a feasibility repair: it may shrink S, but must not enlarge it,
  # since enlarging would change the optimized separation matrix.
  upper = 1                      # maximum permitted gamma
  V_rate = (nsko + 1)/nsko
  for (j in seq_len(num_iter)) {
    gamma = (lower + upper) / 2
    V = V_rate * Sigma - gamma * S
    ok = FALSE
    # Try Cholesky to check PSD w/ margin tol
    try({
      chol(V - tol * diag(nrow(V)))
      ok = TRUE
    }, silent = TRUE)
    if (ok) {
      lower = gamma
    } else {
      upper = gamma
    }
  }
  S_scaled = lower * S * 1   # conservative scale
  list(S = S_scaled, gamma = lower)
}
mvr_loss = function(Sigma, S, smoothing = 0) {
  eigs_S = diag(S)
  eigs_diff = eigen(2 * Sigma - S, symmetric = TRUE, only.values = TRUE)$values

  if (min(eigs_S) < 0 || min(eigs_diff) < 0) {
    return(Inf)
  }

  sum(1 / (eigs_diff + smoothing)) + sum(1 / (eigs_S + smoothing))
}
maxent_loss = function(Sigma, S, I, nsko = 1, smoothing = 1e-05) {
  eigs_S = diag(S[I, I, drop = FALSE])
  V_rate = (nsko + 1)/nsko
  eigs_diff = eigen(V_rate * Sigma - S, symmetric = TRUE, only.values = TRUE)$values
  # smoothing to take care of numerical instability
  det_invG = sum(log(1/(eigs_diff + smoothing))) +
    sum(log(1/(eigs_S + smoothing)))
  det_invG
}
