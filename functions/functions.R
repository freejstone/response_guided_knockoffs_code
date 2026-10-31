#libraries
library(pracma)
library(foreach)
library(doParallel)
library(glmnet)
library(doRNG)
library(EnvStats)
library(stepdownfdp)


###########
#functions
###########
#chol_rank_one_update: Rank-one update/downdate of an upper Cholesky factor.
#solve_mvr_quadratic: Coordinate update used by the MVR optimiser.
#solve_mvr_cd: Ungrouped MVR coordinate-descent solver. (based off knockpy.mrc.solve_mvr)
#solve_maxent_sdp_cd: Solver for entropy knockoffs. (based off knockpy.mrc._solve_maxent_sdp_cd)
#TDC_flex_c: Performs SSS+ or more generally ASS. (based off https://arxiv.org/pdf/1911.09442)
#create_candidate_batch_G: Returns the gram matrix G for a given vector s0. 
#optimal_batch_G: Determines the optimal G using s_equi.
#get_winning_scores: Get_winning_scores using mirandom (based off https://arxiv.org/pdf/1911.09442)
#get_z_scores_batch: Get scores for each feature and knockoff in the same batch.
#create_batch_ko_s0: Get knockoffs for a given s0.
#create_fixed_equi_sdp: Get full-design fixed-X equicorrelated or SDP knockoffs.
#create_fixed_maxent: Get full-design fixed-X maximum-entropy knockoffs.
#create_fixed_mvr: Get full-design fixed-X minimum-variance reconstructability knockoffs.
#augment_fixed_x: Pad p < n < 2p fixed-X problems using the knockoff package convention.
#response_guided_kf: Runs response guided knockoffs
#signed_kf: Runs signed knockoffs
#response_guided_kf_FDX: Runs response guided knockoffs with FDX control
#js_kfwer: Runs JS procedure
#response_guided_kf_kFWER: Runs response guided knockoffs with kFWER control

chol_rank_one_update = function(R, x, add = TRUE) {
  R = as.matrix(R)
  x = as.numeric(x)
  p = length(x)

  if (!all(dim(R) == c(p, p))) {
    stop('R must be a square matrix with dimension length(x)')
  }

  for (k in seq_len(p)) {
    old_diag = R[k, k]
    radicand = old_diag^2 + if (add) x[k]^2 else -x[k]^2
    roundoff_tol = 100 * .Machine$double.eps * max(1, old_diag^2, x[k]^2)
    if (!is.finite(radicand) || radicand < -roundoff_tol) {
      stop('Cholesky downdate is not positive definite')
    }
    radicand = max(radicand, 0)
    new_diag = sqrt(radicand)
    if (old_diag <= 0 || new_diag <= 0) {
      stop('Cholesky update produced a non-positive diagonal')
    }

    c_k = new_diag / old_diag
    s_k = x[k] / old_diag
    R[k, k] = new_diag

    if (k < p) {
      tail = (k + 1):p
      if (add) {
        R[k, tail] = (R[k, tail] + s_k * x[tail]) / c_k
      } else {
        R[k, tail] = (R[k, tail] - s_k * x[tail]) / c_k
      }
      x[tail] = c_k * x[tail] - s_k * R[k, tail]
    }
  }

  R
}

solve_mvr_quadratic = function(cn, cd, sj, smoothing = 0) {
  # Minimise 1/(sj + delta) - delta*cn/(1 - delta*cd), as in knockpy.
  coef2 = -cn - cd^2
  coef1 = 2 * (-cn * (sj + smoothing) + cd)
  coef0 = -cn * (sj + smoothing)^2 - 1
  roots = polyroot(c(coef0, coef1, coef2))

  real_tol = 1e-8 * pmax(1, abs(Re(roots)))
  options = Re(roots[is.finite(Re(roots)) & is.finite(Im(roots)) & abs(Im(roots)) <= real_tol])
  upper_bound = 1 / cd
  lower_bound = -sj
  options = options[options < upper_bound & options > lower_bound]

  if (length(options) == 0) {
    stop(sprintf(
      'All MVR coordinate roots were infeasible or complex (cn=%.6g, cd=%.6g, sj=%.6g)',
      cn, cd, sj
    ))
  }

  losses = 1 / (sj + options) - (options * cn) / (1 - options * cd)
  options[which.min(losses)]
}

solve_mvr_cd = function(Sigma, tol = 1e-8, verbose = TRUE,
                        num_iter = 50, converge_tol = 1e-3,
                        smoothing = 0) {
  V = as.matrix(Sigma)
  if (nrow(V) != ncol(V) || any(!is.finite(V))) {
    stop('Sigma must be a finite square matrix')
  }
  V = (V + t(V)) / 2
  p = ncol(V)
  min_eig = calc_mineig(V)
  if (min_eig <= 0) {
    stop('Sigma must be positive definite for the MVR solver')
  }

  S = diag(rep(min_eig, p), nrow = p, ncol = p)
  R = chol(2 * V - S + smoothing * diag(p))
  loss = Inf
  decayed_improvement = 10
  time0 = proc.time()[[3]]

  for (i in seq_len(num_iter)) {
    for (j in sample.int(p)) {
      ej = numeric(p)
      ej[j] = 1
      vd = forwardsolve(t(R), ej)
      cd = sum(vd^2)
      vn = backsolve(R, vd)
      cn = -sum(vn^2)
      delta = solve_mvr_quadratic(
        cn = cn,
        cd = cd,
        sj = S[j, j],
        smoothing = smoothing
      )

      S[j, j] = S[j, j] + delta
      update_vector = numeric(p)
      update_vector[j] = sqrt(abs(delta))
      R_updated = tryCatch(
        chol_rank_one_update(R, update_vector, add = delta < 0),
        error = function(e) NULL
      )
      if (is.null(R_updated) || any(!is.finite(R_updated))) {
        R = chol(2 * V - S + smoothing * diag(p))
      } else {
        R = R_updated
      }
    }

    # Re-factor once per sweep to prevent accumulated rank-update roundoff.
    R = chol(2 * V - S + smoothing * diag(p))
    previous_loss = loss
    loss = mvr_loss(V, S, smoothing = smoothing)
    if (i > 1) {
      loss_diff = previous_loss - loss
      decayed_improvement = decayed_improvement / 10 + 0.9 * loss_diff
    }
    if (verbose) {
      cat(sprintf(
        'Iter %d | Δt = %.2fs | loss = %.5f | Δ* = %.3e\n',
        i, proc.time()[[3]] - time0, loss, decayed_improvement
      ))
    }
    if (decayed_improvement < converge_tol) {
      if (verbose) cat('Converged – stopping.\n')
      break
    }
  }

  S = shift_until_PSD(S, tol = tol)
  scale_until_PSD(V, S, tol = tol, num_iter = 10)$S
}

solve_maxent_sdp_cd = function(Sigma, I, nsko = 1, tol = 1e-10, verbose = TRUE,
                               num_iter = 50, converge_tol = 1e-4,
                               lambd = 0.5) {
  V   = Sigma
  p   = ncol(V)
  I   = as.integer(I)
  loss = Inf
  
  mineig = calc_mineig(V)
  S_vec  = rep(0, p)
  S_vec[I] = mineig
  S = diag(S_vec, p)
  
  V_rate = (nsko + 1)/nsko
  
  # Initial Cholesky of (2V − S)
  L = chol(V_rate * V - S)
  lambd = min(V_rate * mineig, lambd)
  
  decayed_improvement = 1
  time0 = proc.time()[[3]]
  
  for (i in seq_len(num_iter)) {
    I = sample(I)   # shuffle coordinate order
    for (j in I) {
      diff = V_rate * V - S
      tildey = V_rate * V[j, ]
      tildey[j] = 0
      x = forwardsolve(t(L), tildey)  # solve L x = tildey (L is upper‑tri in R’s chol)
      zeta = diff[j, j]
      x22  = sum(x^2)
      qinv = zeta * x22 / (zeta + x22)
      sj_star = (V_rate * V[j, j] - qinv) / 2
      # Update S and re‑factorise – slower but stable
      S[j, j] = sj_star
      L = chol(V_rate * V - S)
    }
    # Convergence metrics
    prev_loss = loss
    loss = maxent_loss(V, S, I, nsko = nsko)
    if (i > 1) {
      loss_diff = prev_loss - loss
      decayed_improvement = decayed_improvement / 10 + 0.9 * loss_diff
    }
    if (verbose) {
      cat(sprintf("Iter %d | Δt = %.2fs | loss = %.5f | Δ* = %.3e\n",
                  i, proc.time()[[3]] - time0, loss, decayed_improvement))
    }
    if (decayed_improvement < converge_tol) {
      if (verbose) cat("Converged – stopping.\n")
      break
    }
  }
  
  # Final feasibility shifts & scaling
  S[I, I] = shift_until_PSD(S[I, I], tol = tol)
  scaled  = scale_until_PSD(V, S, nsko = nsko, tol = tol, num_iter = 50)
  scaled$S
}
TDC_flex_c = function(KF, obs, score = NULL, BC1 = 1, c = 1/2, lambda = 1/2) {
  if (!is.null(score)) {
    indxs = order(score, decreasing = TRUE)
    KF = KF[indxs]
    obs = obs[indxs]
  }
  nTD = cumsum(obs)
  nDD = cumsum(KF)
  fdps = (pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals = rev(cummin(rev(fdps)))
  if (!is.null(score)) {
    qvals[indxs] = qvals
  }
  return(qvals)
}
create_candidate_batch_G = function(I, Sig, nsko, s0){
  p = ncol(Sig)
  lenI = length(I)
  lenG = nsko*lenI + p
  G = matrix(0, lenG, lenG)
  P = 1:p
  if (length(s0) == 1){
    diag_s = diag(rep(s0, p), nrow = p, ncol = p)
    Sig0 = Sig - diag_s
  } else if (length(s0) > 1){
    diag_s = diag(c(s0))
    I_IP = matrix(0, nrow = lenI, ncol = p)
    I_IP[(I - 1)*lenI + (1:lenI)] = 1
    I_PI = t(I_IP)
    Sig0 = Sig - I_PI%*%diag_s%*%I_IP
  }
  SigII = Sig[I,I]
  Sig0II = Sig0[I,I]
  Sig0IP = Sig0[I,P]
  Sig0PI = Sig0[P,I]
  
  G[P, P] = Sig
  for (i in 1:nsko){
    G[P, (p+(i-1)*lenI+1) : (p+i*lenI)] = Sig0PI
    G[(p+(i-1)*lenI+1) : (p+i*lenI), P] = Sig0IP
  }
  for (i in 1:nsko){
    G[(p+(i-1)*lenI+1) : (p+i*lenI), (p+(i-1)*lenI+1) : (p+i*lenI)] = SigII
    for (j in 1:nsko){
      if (j != i){
        G[(p+(i-1)*lenI+1) : (p+i*lenI) , (p+(j-1)*lenI+1) : (p+j*lenI)] = Sig0II
      }
    }
  }
  return(G)
}
optimal_batch_G = function(I, Sig, nsko, analytic = TRUE, s0_only = FALSE){
  minEig = function(s){
    mineig = min(eigen(create_candidate_batch_G(I, Sig, nsko, s))$values)
    if (abs(mineig) < 1e-13){
      mineig = 0
    }
    return(mineig)
  }
  sigeig = eigen(Sig)$values
  s0 = min(sigeig)
  tol = 1e-5
  if (any(sigeig <= tol*max(sigeig))){
    warning('Model is not identifiable, but proceeding with knockoffs')
  }
  # analytical result says that s_equi for NO batching is sigeig*(nsko + 1)/nsko
  # so, batching has to be at least as good as this
  if (minEig(s0*(nsko + 1)/nsko) == 0){
    s_equi = s0*(nsko + 1)/nsko
  } else {
    if (analytic) {
      P = setdiff(1:nrow(Sig), I)
      s_equi = ((nsko + 1)/nsko)*min(eigen(Sig[I, I] - Sig[I, P] %*% pracma::inv(Sig[P, P]) %*% Sig[P, I])$values)
    } else {
      s_equi = as.numeric(fzero(minEig, c(s0*(nsko + 1)/nsko, 2))$x)
    }
  }
  if (s0_only) {
    G = NULL
  } else {
    G = create_candidate_batch_G(I, Sig, nsko, min(1, s_equi))
  }
  return(list(G = G, s_equi = s_equi))
}
get_winning_scores = function(obs_scores, decoys, c, lambda){
  obs_scores = as.matrix(obs_scores)
  decoys = as.matrix(decoys)
  m = dim(decoys)[1]
  n_decoys = dim(decoys)[2]
  mixed_scores = -.Machine$double.xmax * ones(m, 1)
  
  sorted_decoys = t(apply(decoys,1,sort))
  n_decoys_LT_obs = as.matrix(rowSums(sweep(decoys, 1, obs_scores) < 0))
  n_decoys_EQ_obs = as.matrix(rowSums(sweep(decoys, 1, obs_scores) == 0))
  n_decoys_GT_obs = n_decoys - (n_decoys_LT_obs + n_decoys_EQ_obs)
  n_decoys_NE_obs = n_decoys - n_decoys_EQ_obs
  p_obs = (n_decoys_NE_obs - n_decoys_LT_obs + 1) / (n_decoys_NE_obs + 1)
  obs_wins = as.matrix(p_obs <= c)
  mixed_scores[which(obs_wins == 1)] = obs_scores[which(obs_wins == 1)]
  
  #	Define the mapping of the losing observation ranks to the winning decoys ranks
  n_decoys_in_draw = floor(c*(n_decoys + 1)) # the number of decoys the losing observations will be mapped to
  n_obs_ranks = n_decoys - ceil((n_decoys + 1) * lambda) + 1	# the number of ranks in which the observation loses
  max_mapped_decoy_rank = zeros(n_obs_ranks, 1) # the highest decoy rank each losing observation is mapped to
  max_mapped_decoy_prob = zeros(n_obs_ranks, 1) # the probability of mapping to that maximal decoy rank
  current_decoy_rank = n_decoys # the decoy rank currently being mapped into, starting from the highest rank
  current_decoy_coverage = 0		# how much of that decoy rank did we cover so far
  uni_decoy_coverage = n_obs_ranks / n_decoys_in_draw	# how much coverage should each mapped-to-decoy get
  for (i in 1 : n_obs_ranks){ # loop defining the map, one observed rank at a time
    max_mapped_decoy_rank[i] = current_decoy_rank		# assign the current highest decoy rank available
    if (current_decoy_coverage + 1 > uni_decoy_coverage){ # mapping the current observation (i) to this decoy overfills its quota
      max_mapped_decoy_prob[i] = uni_decoy_coverage - current_decoy_coverage	# so the probability is (quota) - (coverage so far)
      remainder_obs_coverage_prob = 1 - max_mapped_decoy_prob[i]		# portion of the current observation we still need to map
      current_decoy_coverage = uni_decoy_coverage;				# current decoy is saturated
    } else { # current observation can be fully mapped to current decoy
      max_mapped_decoy_prob[i] = 1					# so the probability is 1
      current_decoy_coverage = 1 + current_decoy_coverage	# update current decoy coverage
      remainder_obs_coverage_prob = 0				# there is no "change", the current observation was fully mapped
    }
    if (current_decoy_coverage >= uni_decoy_coverage - 1e-12){ # is current decoy saturated (up to a roundoff error)?
      # floor(...) below is the number of additional decoy ranks the current observation can saturate
      current_decoy_rank = current_decoy_rank - floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) - 1
      current_decoy_coverage = remainder_obs_coverage_prob - floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) * uni_decoy_coverage	# and the remainder of the observation is allocated to the new decoy
    }
  }
  if (n_decoys_in_draw > 0 && (current_decoy_coverage > 1e-12 || n_decoys - current_decoy_rank != n_decoys_in_draw)){	# sanity check - increase 1e-12 for larger m
    stop('Warning, the mapping doesnt work...')
  }
  # next, apply the mapping to every observation that lost
  decoy_wins = (n_decoys_in_draw > 0) & (n_decoys_GT_obs > 0) & !obs_wins
  counted_decoy_wins = decoy_wins & (p_obs > lambda)
  obs_ranks = n_decoys_LT_obs[counted_decoy_wins] + floor((n_decoys_EQ_obs[counted_decoy_wins]+1)*rand(sum(counted_decoy_wins), 1))	# randomly broken ties
  obs_ranks = pmin(obs_ranks, n_obs_ranks-1)		# because ties are randomly broken the initial rank can exceed n_obs_ranks-1
  # 	if c <= 0.5 mapped_obs_ranks = max_mapped_decoy_rank(1+obs_ranks) - (rand(sum(decoy_wins), 1) > max_mapped_decoy_prob(1+obs_ranks));
  rands = rand(sum(counted_decoy_wins), 1)
  mapped_obs_ranks = max_mapped_decoy_rank[1+obs_ranks] - ceil((rands - max_mapped_decoy_prob[1+obs_ranks]) / uni_decoy_coverage - 1e-12)*(rands > max_mapped_decoy_prob[1+obs_ranks])		# this is a bit tricky but it works (look at the notes to see why)
  mixed_scores[which(counted_decoy_wins == T)] = sorted_decoys[(c(mapped_obs_ranks) - 1)*nrow(sorted_decoys) + which(counted_decoy_wins == T)] # the selected scores
  
  return(list(decoy_wins = decoy_wins, counted_decoy_wins = counted_decoy_wins, obs_wins = obs_wins, mixed_scores = mixed_scores, p_obs = p_obs))
}
get_z_scores_batch = function(X, X_k, y, I, nsko) {
  p = ncol(X)
  n = nrow(X)
  X_all = cbind(X, X_k)
  swap = sample(1:ncol(X_all))
  Zs = knockoff:::lasso_max_lambda(X_all[, swap], y, intercept = F, standardize = F, nlambda = 2000)
  Zs[swap] = Zs
  Z = Zs[I]
  Zk = matrix(0, nrow = length(I), ncol = nsko)
  for (i in 1:nsko) {
    Zk[, i] = Zs[ ((p + 1) + (i - 1)*length(I)): ((p + length(I)) + (i - 1)*length(I))]
  }
  return(list(Z = Z, Zk = Zk))
}
create_batch_ko_s0 = function(X, I, nsko, s0 = NULL){
  n = nrow(X)
  p = ncol(X)
  lenI = length(I)
  
  thinX = cbind(X, pracma::randn(n, nsko*lenI))
  QRthinX = qr(thinX)
  Qall = qr.Q(QRthinX)
  Rall = qr.R(QRthinX)[, QRthinX$pivot]
  diagRall = diag(Rall)
  
  Sig = t(X)%*%X
  if (!any(is.null(s0))) {
    G = create_candidate_batch_G(I, Sig, nsko, s0)
  } else {
    Sig = t(X)%*%X
    Results = optimal_batch_G(I, Sig, nsko)
    G = Results[[1]] #batch specific Gram
    s0 = Results[[2]] #corresponding s0
  }
  dimG = nsko*lenI + p #p + d|I|
  Decom = eigen(G, symmetric = TRUE)
  X0 = Decom$vectors %*% diag(sqrt(pmax(Decom$values, 0))) %*% t(Decom$vectors) #construction of X0 through SVD
  qr_X0 = qr(X0)
  Q0 = qr.Q(qr_X0)
  R0 = qr.R(qr_X0)
  R0[, qr_X0$pivot] = R0
  
  QI = Qall
  diagRI = diagRall[1:p]
  
  #makes R0 agree with R so that first p of X1 is X!
  signs = cbind(
    matrix(sign(diagRI / diag(R0[1:p, 1:p, drop = FALSE])), nrow = 1),
    pracma::ones(1, nsko * lenI)
  )
  Qs = QI %*% diag(c(signs))
  X1 = Qs %*% R0
  
  X_ko = X1[, (p+1):dimG]
  Xall = cbind(X, X_ko)
  if (any( abs(t(Xall)%*%Xall - G) > 1e-7)){
    print(max(abs(t(Xall)%*%Xall - G)))
    stop('Xall constructed is different to G (up to 1e-7)')
  }
  return(list(X_ko = X_ko, s0 = s0))
}
create_fixed_equi_sdp = function(X, type = c('equi', 'sdp'), s0 = NULL,
                                 s0_only = FALSE, verbose = FALSE) {
  type = match.arg(type)
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)

  if (n < 2*p) {
    stop(sprintf(
      'Full-design fixed-X %s knockoffs require n >= 2p; got n = %d and p = %d',
      toupper(type), n, p
    ))
  }
  if (any(!is.finite(X))) {
    stop('X must contain only finite values')
  }

  column_norms = sqrt(colSums(X^2))
  if (any(column_norms <= sqrt(.Machine$double.eps))) {
    stop('X contains a zero or numerically zero column')
  }
  X = sweep(X, 2, column_norms, '/')
  if (qr(X)$rank < p) {
    stop(sprintf('X must have full column rank for full-design fixed-X %s knockoffs', toupper(type)))
  }

  Sigma = t(X)%*%X
  if (is.null(s0)) {
    s0 = switch(
      type,
      equi = {
        singular_values = svd(X, nu = 0, nv = 0)$d
        if (any(singular_values <= 1e-5 * max(singular_values))) {
          stop('Data matrix is rank deficient. Equicorrelated knockoffs will have no power.')
        }
        rep(min(2 * min(singular_values)^2, 1), p)
      },
      sdp = {
        solved_s0 = knockoff::create.solve_sdp(Sigma, verbose = verbose)
        # Match knockoff::create.fixed(method = 'sdp'), which removes tiny entries
        # before constructing the randomized complement.
        solved_s0[solved_s0 <= 1e-5] = 0
        solved_s0
      }
    )
  }
  if (length(s0) != p || any(!is.finite(s0)) || any(s0 < -1e-10) || any(s0 > 1 + 1e-7)) {
    stop(sprintf('The %s solver returned an invalid separation vector', toupper(type)))
  }
  s0 = pmax(s0, 0)
  S = diag(s0, nrow = p, ncol = p)
  if (min(eigen(2*Sigma - S, symmetric = TRUE, only.values = TRUE)$values) < -1e-7) {
    stop(sprintf('The %s separation vector is not feasible', toupper(type)))
  }

  if (s0_only) {
    return(list(X = X, s0 = s0))
  }

  # Mirror knockoff::create.fixed after its separation-vector calculation. This
  # preserves its randomized construction and RNG consumption while allowing the
  # equicorrelated/SDP optimization itself to be cached.
  X_svd = knockoff:::decompose(X, randomize = TRUE)
  if (type == 'equi') {
    if (max(abs(s0 - s0[1])) > 1e-10) {
      stop('The equicorrelated separation vector must be constant')
    }
    s = s0[1]
    s_diff = pmax(0, 2*s - (s/X_svd$d)^2)
    Xk = (
      X_svd$u %*% diag(X_svd$d - s/X_svd$d, nrow = p, ncol = p) +
        X_svd$u_perp %*% diag(sqrt(s_diff), nrow = p, ncol = p)
    ) %*% t(X_svd$v)
  } else {
    d_inv = 1/X_svd$d
    d_inv[X_svd$d <= 1e-5 * max(X_svd$d)] = 0
    G_inv = (
      X_svd$v %*% diag(d_inv^2, nrow = p, ncol = p)
    ) %*% t(X_svd$v)
    S = diag(s0, nrow = p, ncol = p)
    C_svd = knockoff:::canonical_svd(2*S - S %*% G_inv %*% S)
    Xk = X - X %*% G_inv %*% S + (
      X_svd$u_perp %*% diag(sqrt(pmax(0, C_svd$d)), nrow = p, ncol = p)
    ) %*% t(C_svd$v)
  }

  return(list(X = X, Xk = Xk, s0 = s0))
}
create_fixed_maxent = function(X, s0 = NULL, s0_only = FALSE, verbose = FALSE) {
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)

  if (n < 2*p) {
    stop(sprintf('Full-design fixed-X ME knockoffs require n >= 2p; got n = %d and p = %d', n, p))
  }
  if (any(!is.finite(X))) {
    stop('X must contain only finite values')
  }

  column_norms = sqrt(colSums(X^2))
  if (any(column_norms <= sqrt(.Machine$double.eps))) {
    stop('X contains a zero or numerically zero column')
  }
  X = sweep(X, 2, column_norms, '/')
  if (qr(X)$rank < p) {
    stop('X must have full column rank for full-design fixed-X ME knockoffs')
  }

  Sigma = t(X)%*%X
  if (is.null(s0)) {
    S = solve_maxent_sdp_cd(Sigma, I = 1:p, verbose = verbose)
    s0 = diag(S)
  }
  if (length(s0) != p || any(!is.finite(s0)) || any(s0 < -1e-10) || any(s0 > 1 + 1e-7)) {
    stop('The maximum-entropy solver returned an invalid separation vector')
  }
  s0 = pmax(s0, 0)
  if (min(eigen(2*Sigma - diag(s0, p), symmetric = TRUE, only.values = TRUE)$values) < -1e-7) {
    stop('The maximum-entropy separation vector is not feasible')
  }

  if (s0_only) {
    return(list(X = X, s0 = s0))
  }

  X_ko = create_batch_ko_s0(X, 1:p, nsko = 1, s0 = s0)
  return(list(X = X, Xk = X_ko$X_ko, s0 = s0))
}
create_fixed_mvr = function(X, s0 = NULL, s0_only = FALSE,
                            verbose = FALSE, tol = NULL) {
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)

  if (n < 2*p) {
    stop(sprintf('Full-design fixed-X MVR knockoffs require n >= 2p; got n = %d and p = %d', n, p))
  }
  if (any(!is.finite(X))) {
    stop('X must contain only finite values')
  }

  column_norms = sqrt(colSums(X^2))
  if (any(column_norms <= sqrt(.Machine$double.eps))) {
    stop('X contains a zero or numerically zero column')
  }
  X = sweep(X, 2, column_norms, '/')
  if (qr(X)$rank < p) {
    stop('X must have full column rank for full-design fixed-X MVR knockoffs')
  }

  Sigma = t(X)%*%X
  if (is.null(tol)) {
    # Match knockpy's fixed-X tolerance, DEFAULT_TOL / n.
    tol = 1e-5 / n
  }
  if (is.null(s0)) {
    S = solve_mvr_cd(Sigma, tol = tol, verbose = verbose)
    s0 = diag(S)
  }
  if (length(s0) != p || any(!is.finite(s0)) || any(s0 < -1e-10)) {
    stop('The MVR solver returned an invalid separation vector')
  }
  s0 = pmax(s0, 0)
  if (min(eigen(2*Sigma - diag(s0, p), symmetric = TRUE, only.values = TRUE)$values) < -1e-7) {
    stop('The MVR separation vector is not feasible')
  }

  if (s0_only) {
    return(list(X = X, s0 = s0))
  }

  X_ko = create_batch_ko_s0(X, 1:p, nsko = 1, s0 = s0)
  return(list(X = X, Xk = X_ko$X_ko, s0 = s0))
}
augment_fixed_x = function(X, y) {
  X = as.matrix(X)
  y = as.numeric(y)
  n = nrow(X)
  p = ncol(X)

  if (n <= p) {
    stop(sprintf('Fixed-X knockoffs require n > p; got n = %d and p = %d', n, p))
  }
  if (any(!is.finite(X))) {
    stop('X must contain only finite values')
  }
  if (length(y) != n || any(!is.finite(y))) {
    stop('y must be a finite vector with one entry per row of X')
  }
  if (n >= 2*p) {
    return(list(X = X, y = y))
  }

  X_qr = qr(X)
  if (X_qr$rank < p) {
    stop('X must have full column rank for fixed-X knockoffs')
  }
  residuals = qr.resid(X_qr, y)
  sigma = sqrt(sum(residuals^2) / (n - p))
  n_extra = 2*p - n

  # This is the same randomized augmentation used by knockoff::create.fixed:
  # append zero design rows and responses drawn using the OLS residual scale.
  X = rbind(X, matrix(0, nrow = n_extra, ncol = p))
  y = c(y, rnorm(n_extra, sd = sigma))
  list(X = X, y = y)
}
response_guided_kf = function(y, y_probe, X, pos, beta = NULL, 
                              d = 1, type = 'entropy', 
                              alphalist = seq(0.3, 0.01, by = -0.01),
                              seed = 1) {
  set.seed(seed)
  n = nrow(X)
  p = ncol(X)
  
  # Extend
  Xe = cbind(X, y_probe, rep(1, n))
  Xe = apply(Xe, 2, function(x) x/sqrt(sum(x^2)))
  X = apply(X, 2, function(x) x/sqrt(sum(x^2)))
  Sigma.e = t(Xe)%*%Xe
  
  # Construct knockoffs
  if (type == 'entropy') {
    C = setdiff(1:(p + 2), pos)
    Sigma.e.reduced = Sigma.e[pos, pos] - Sigma.e[pos, C] %*% inv(Sigma.e[C, C]) %*% Sigma.e[C, pos]
    s0 = diag(solve_maxent_sdp_cd(Sigma.e.reduced, I = 1:ncol(Sigma.e.reduced)))
    s0 = s0[s0 > 0]
  } else if (type == 'equi') {
    s0.res = optimal_batch_G(pos, Sigma.e, nsko = 1, s0_only = TRUE)
    s0 = rep(s0.res$s_equi, length(pos))
  } 
  X_ko = create_batch_ko_s0(Xe, pos, nsko = 1, s0 = s0)
  
  # Score and signs
  scores = get_z_scores_batch(X, X_ko$X_ko, y, pos, 1)
  W = pmax(scores$Z, scores$Zk)
  W[scores$Z < scores$Zk] = - W[scores$Z < scores$Zk]
  W[scores$Z == scores$Zk] = 0
  signs = sign(t(X[, pos] - X_ko$X_ko) %*%y)*sign(W)
  
  # FDR analysis
  qvals = TDC_flex_c(W[signs*d > 0] < 0, W[signs*d > 0] > 0, score = abs(W[signs*d > 0]), c = 1/2, lambda = 1/2)
  discs = vector(mode = "list", length = length(alphalist))
  for (i in 1:length(alphalist)) {
    disc = (qvals <= alphalist[i]) & (W[signs*d > 0] > 0)
    discs[[i]] = disc
  }
  if (is.null(beta)) return(list(discs = discs))
  power = rep(0, length(alphalist))
  fdp = rep(0, length(alphalist))
  for (i in 1:length(alphalist)) {
    disc = (qvals <= alphalist[i]) & (W[signs*d > 0] > 0)
    power[i] = sum(disc & (beta[pos][signs*d > 0]*d > 0))/max(sum(beta*d > 0), 1)
    fdp[i] = sum(disc & (beta[pos][signs*d > 0]*d <= 0))/max(sum(disc), 1)
  }
  return(list(discs = discs, power = power, fdp = fdp))
}
signed_kf = function(y, X, beta, d = 1, type = 'equi', alphalist = seq(0.3, 0.01, by = -0.01), seed = 1, s0 = NULL) {
  set.seed(seed)
  n = nrow(X)
  p = ncol(X)

  if (type %in% c('entropy', 'maxent', 'me')) {
    augmented = augment_fixed_x(X, y)
    X = augmented$X
    y = augmented$y
    X_ko = create_fixed_maxent(X, s0 = s0, verbose = FALSE)
    X = X_ko$X
    Xk = X_ko$Xk
  } else if (type == 'mvr') {
    augmented = augment_fixed_x(X, y)
    X = augmented$X
    y = augmented$y
    X_ko = create_fixed_mvr(X, s0 = s0, verbose = FALSE)
    X = X_ko$X
    Xk = X_ko$Xk
  } else if (type %in% c('equi', 'sdp') && !is.null(s0)) {
    augmented = augment_fixed_x(X, y)
    X = augmented$X
    y = augmented$y
    X_ko = create_fixed_equi_sdp(X, type = type, s0 = s0, verbose = FALSE)
    X = X_ko$X
    Xk = X_ko$Xk
  } else {
    X_ko = knockoff::create.fixed(X, method = type, y = y, randomize = TRUE)
    X = X_ko$X
    Xk = X_ko$Xk
    y = X_ko$y
  }

  scores = get_z_scores_batch(X, Xk, y, 1:p, 1)
  W = pmax(scores$Z, scores$Zk)
  W[scores$Z < scores$Zk] = -W[scores$Z < scores$Zk]
  W[scores$Z == scores$Zk] = 0
  W.pos.KF = sum(W > 0 & beta != 0)
  signs = sign(t(X - Xk) %*%y)*sign(W)
  
  # FDR analysis
  qvals = TDC_flex_c(W[signs*d > 0] < 0, W[signs*d > 0] > 0, score = abs(W[signs*d > 0]), c = 1/2, lambda = 1/2)
  discs = vector(mode = "list", length = length(alphalist))
  for (i in 1:length(alphalist)) {
    disc = (qvals <= alphalist[i]) & (W[signs*d > 0] > 0)
    discs[[i]] = disc
  }
  if (is.null(beta)) return(list(discs = discs))
  power = rep(0, length(alphalist))
  fdp = rep(0, length(alphalist))
  for (i in 1:length(alphalist)) {
    disc = discs[[i]]
    power[i] = sum(disc & (beta[signs*d > 0]*d > 0))/max(sum(beta*d > 0), 1)
    fdp[i] = sum(disc & (beta[signs*d > 0]*d <= 0))/max(sum(disc), 1)
  }
  return(list(discs = discs, power = power, fdp = fdp))
}
response_guided_kf_FDX = function(y, y_probe, X, pos, beta = NULL, 
                              d = 1, type = 'entropy', 
                              alphalist = seq(0.3, 0.01, by = -0.01),
                              gammalist = c(0.5, 0.1, 0.05),
                              seed = 1) {
  set.seed(seed)
  n = nrow(X)
  p = ncol(X)
  
  # Extend
  Xe = cbind(X, y_probe, rep(1, n))
  Xe = apply(Xe, 2, function(x) x/sqrt(sum(x^2)))
  X = apply(X, 2, function(x) x/sqrt(sum(x^2)))
  Sigma.e = t(Xe)%*%Xe
  
  # Construct knockoffs
  if (type == 'both') {
    pos_sign = which(t(X) %*% y > 0)
    
    # Get entropy
    C = setdiff(1:p, pos_sign)
    Sigma.e.reduced = Sigma.e[pos_sign, pos_sign] - Sigma.e[pos_sign, C] %*% inv(Sigma.e[C, C]) %*% Sigma.e[C, pos_sign]
    s0.entropy = diag(solve_maxent_sdp_cd(Sigma.e.reduced, I = 1:ncol(Sigma.e.reduced)))
    s0.entropy = s0.entropy[s0.entropy > 0]
    
    # Get equi
    s0.res = optimal_batch_G(pos_sign, Sigma.e[1:p, 1:p], nsko = 1, s0_only = TRUE)
    s0.equi = rep(s0.res$s_equi, length(pos_sign))
    
    # Try each of them and decide on the better one
    X_ko.entropy = create_batch_ko_s0(X, pos_sign, nsko = 1, s0 = s0.entropy)
    X_ko.equi = create_batch_ko_s0(X, pos_sign, nsko = 1, s0 = s0.equi)
    
    # Score 
    scores = get_z_scores_batch(X, X_ko.entropy$X_ko, y, pos_sign, 1)
    W = pmax(scores$Z, scores$Zk)
    W[scores$Z < scores$Zk] = - W[scores$Z < scores$Zk]
    W[scores$Z == scores$Zk] = 0
    signs = sign(t(X[, pos_sign] - X_ko.entropy$X_ko) %*%y)*sign(W)
    
    # FDR analysis
    qvals = TDC_flex_c(W[signs*d > 0] < 0, W[signs*d > 0] > 0, score = abs(W[signs*d > 0]), c = 1/2, lambda = 1/2, BC1 = 0)
    disc.entropy = sum((qvals <= 0.5) & (W[signs*d > 0] > 0))
    
    
    # Score 
    scores = get_z_scores_batch(X, X_ko.equi$X_ko, y, pos_sign, 1)
    W = pmax(scores$Z, scores$Zk)
    W[scores$Z < scores$Zk] = - W[scores$Z < scores$Zk]
    W[scores$Z == scores$Zk] = 0
    signs = sign(t(X[, pos_sign] - X_ko.equi$X_ko) %*%y)*sign(W)
    
    # FDR analysis
    qvals = TDC_flex_c(W[signs*d > 0] < 0, W[signs*d > 0] > 0, score = abs(W[signs*d > 0]), c = 1/2, lambda = 1/2, BC1 = 0)
    disc.equi = sum((qvals <= 0.2) & (W[signs*d > 0] > 0))
    
    print(disc.equi)
    print(disc.entropy)
    if (disc.equi >= disc.entropy) {
      type = 'equi'
    } else {
      type = 'entropy'
    }
    print(paste('choose:', type))
  }
  
  if (type == 'entropy') {
    C = setdiff(1:(p + 2), pos)
    Sigma.e.reduced = Sigma.e[pos, pos] - Sigma.e[pos, C] %*% inv(Sigma.e[C, C]) %*% Sigma.e[C, pos]
    s0 = diag(solve_maxent_sdp_cd(Sigma.e.reduced, I = 1:ncol(Sigma.e.reduced)))
    s0 = s0[s0 > 0]
  } else if (type == 'equi') {
    s0.res = optimal_batch_G(pos, Sigma.e, nsko = 1, s0_only = TRUE)
    s0 = rep(s0.res$s_equi, length(pos))
  } 
  X_ko = create_batch_ko_s0(Xe, pos, nsko = 1, s0 = s0)
  
  # Score and signs
  scores = get_z_scores_batch(X, X_ko$X_ko, y, pos, 1)
  W = pmax(scores$Z, scores$Zk)
  W[scores$Z < scores$Zk] = - W[scores$Z < scores$Zk]
  W[scores$Z == scores$Zk] = 0
  signs = sign(t(X[, pos] - X_ko$X_ko) %*%y)*sign(W)
  
  # FDX analysis
  keep = which(signs*d > 0)
  W_d = W[keep]
  scores_and_labels = cbind(abs(W_d), sign(W_d))
  discs = list()
  count = 1
  for (alpha in alphalist) {
    for (gamma in gammalist) {
      if (length(keep) == 0) {
        res = list(discoveries = numeric(0), discoveries_ind = integer(0))
      } else {
        res = fdp_sd(scores_and_labels, alpha, gamma)
      }
      # Map back to original feature indices
      res$discoveries_ind = pos[keep][res$discoveries_ind]
      res$alpha = alpha
      res$gamma = gamma
      if (!is.null(beta)) {
        n_disc = length(res$discoveries_ind)
        n_true = sum(beta[res$discoveries_ind]*d > 0)
        res$power = n_true/max(sum(beta*d > 0), 1)
        res$fdp = (n_disc - n_true)/max(n_disc, 1)
      }
      discs[[count]] = res
      count = count + 1
    }
  }
  return(list(discs = discs))
}
js_kfwer = function(W, k, alpha, max_w = 1e6) {
  # Check
  if (!is.numeric(W) || any(!is.finite(W))) stop("'W' must be a finite numeric vector.")
  if (length(W) < 1) stop("'W' must have length >= 1.")
  if (length(k) != 1 || k < 1 || k != as.integer(k)) stop("'k' must be a positive integer.")
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be in (0, 1).")
  }
  
  m = length(W)
  
  # 
  # The series sum_{i=k}^\infty 2^{-(i+w)} * choose(i+w-1, i)
  # is exactly P(X >= k) for X ~ NegBin(size=w, prob=1/2) counting failures.
  tail_prob = function(w) {
    stats::pnbinom(q = k - 1, size = w, prob = 0.5, lower.tail = FALSE)
  }
  
  v = 0L
  w = 1L
  while (w <= max_w) {
    if (tail_prob(w) <= alpha) {
      v = w
      w = w + 1L
    } else {
      break
    }
  }
  if (w > max_w) stop("Reached max_w without failing the inequality; increase 'max_w' if needed.")
  
  # If v = 0, Step 2 gives tau = +Inf (effectively), hence no discoveries.
  if (v == 0L) {
    return(list(S_hat = integer(0), tau = Inf, v = v))
  }
  
  # tau = sup{ t > 0 : #{j : W_j <= -t} = v }.
  # With continuous W (no ties), this corresponds to tau = -w_(v),
  # where w_(v) is the v-th smallest entry of W.
  #
  # If there are fewer than v strictly-negative W values, the set is empty; use tau = 0.
  n_neg = sum(W < 0)
  if (v > n_neg) {
    tau = 0
  } else {
    w_sorted = sort(W, decreasing = FALSE)
    tau = -w_sorted[v]
    if (!is.finite(tau) || tau < 0) tau = 0
  }
  
  S_hat = which(W >= tau)
  
  list(S_hat = S_hat, tau = tau, v = v)
}
response_guided_kf_kFWER = function(y, y_probe, X, pos, beta = NULL, 
                                  d = 1, type = 'entropy', 
                                  alphalist = seq(0.3, 0.01, by = -0.01),
                                  klist = c(10, 5, 1),
                                  seed = 1) {
  set.seed(seed)
  n = nrow(X)
  p = ncol(X)
  
  # Extend
  Xe = cbind(X, y_probe, rep(1, n))
  Xe = apply(Xe, 2, function(x) x/sqrt(sum(x^2)))
  X = apply(X, 2, function(x) x/sqrt(sum(x^2)))
  Sigma.e = t(Xe)%*%Xe
  
  # Construct knockoffs
  if (type == 'both') {
    pos_sign = which(t(X) %*% y > 0)
    
    # Get entropy
    C = setdiff(1:p, pos_sign)
    Sigma.e.reduced = Sigma.e[pos_sign, pos_sign] - Sigma.e[pos_sign, C] %*% inv(Sigma.e[C, C]) %*% Sigma.e[C, pos_sign]
    s0.entropy = diag(solve_maxent_sdp_cd(Sigma.e.reduced, I = 1:ncol(Sigma.e.reduced)))
    s0.entropy = s0.entropy[s0.entropy > 0]
    
    # Get equi
    s0.res = optimal_batch_G(pos_sign, Sigma.e[1:p, 1:p], nsko = 1, s0_only = TRUE)
    s0.equi = rep(s0.res$s_equi, length(pos_sign))
    
    # Try each of them and decide on the better one
    X_ko.entropy = create_batch_ko_s0(X, pos_sign, nsko = 1, s0 = s0.entropy)
    X_ko.equi = create_batch_ko_s0(X, pos_sign, nsko = 1, s0 = s0.equi)
    
    # Score 
    scores = get_z_scores_batch(X, X_ko.entropy$X_ko, y, pos_sign, 1)
    W = pmax(scores$Z, scores$Zk)
    W[scores$Z < scores$Zk] = - W[scores$Z < scores$Zk]
    W[scores$Z == scores$Zk] = 0
    signs = sign(t(X[, pos_sign] - X_ko.entropy$X_ko) %*%y)*sign(W)
    
    # FDR analysis
    qvals = TDC_flex_c(W[signs*d > 0] < 0, W[signs*d > 0] > 0, score = abs(W[signs*d > 0]), c = 1/2, lambda = 1/2, BC1 = 0)
    disc.entropy = sum((qvals <= 0.5) & (W[signs*d > 0] > 0))
    
    
    # Score 
    scores = get_z_scores_batch(X, X_ko.equi$X_ko, y, pos_sign, 1)
    W = pmax(scores$Z, scores$Zk)
    W[scores$Z < scores$Zk] = - W[scores$Z < scores$Zk]
    W[scores$Z == scores$Zk] = 0
    signs = sign(t(X[, pos_sign] - X_ko.equi$X_ko) %*%y)*sign(W)
    
    # FDR analysis
    qvals = TDC_flex_c(W[signs*d > 0] < 0, W[signs*d > 0] > 0, score = abs(W[signs*d > 0]), c = 1/2, lambda = 1/2, BC1 = 0)
    disc.equi = sum((qvals <= 0.2) & (W[signs*d > 0] > 0))
    
    print(disc.equi)
    print(disc.entropy)
    if (disc.equi >= disc.entropy) {
      type = 'equi'
    } else {
      type = 'entropy'
    }
    print(paste('choose:', type))
  }
  
  if (type == 'entropy') {
    C = setdiff(1:(p + 2), pos)
    Sigma.e.reduced = Sigma.e[pos, pos] - Sigma.e[pos, C] %*% inv(Sigma.e[C, C]) %*% Sigma.e[C, pos]
    s0 = diag(solve_maxent_sdp_cd(Sigma.e.reduced, I = 1:ncol(Sigma.e.reduced)))
    s0 = s0[s0 > 0]
  } else if (type == 'equi') {
    s0.res = optimal_batch_G(pos, Sigma.e, nsko = 1, s0_only = TRUE)
    s0 = rep(s0.res$s_equi, length(pos))
  } 
  X_ko = create_batch_ko_s0(Xe, pos, nsko = 1, s0 = s0)
  
  # Score and signs
  scores = get_z_scores_batch(X, X_ko$X_ko, y, pos, 1)
  W = pmax(scores$Z, scores$Zk)
  W[scores$Z < scores$Zk] = - W[scores$Z < scores$Zk]
  W[scores$Z == scores$Zk] = 0
  signs = sign(t(X[, pos] - X_ko$X_ko) %*%y)*sign(W)
  
  # k-FWER analysis
  keep = which(signs*d > 0)
  W_d = W[keep]
  discs = list()
  count = 1
  for (alpha in alphalist) {
    for (k in klist) {
      if (length(keep) == 0) {
        res = list(S_hat = integer(0), tau = Inf, v = 0L)
      } else {
        res = js_kfwer(W_d, k = k, alpha = alpha)
      }
      # Map back to original feature indices
      res$S_hat = pos[keep][res$S_hat]
      res$alpha = alpha
      res$k = k
      if (!is.null(beta)) {
        n_disc = length(res$S_hat)
        n_true = sum(beta[res$S_hat]*d > 0)
        res$power = n_true/max(sum(beta*d > 0), 1)
        res$fdp = (n_disc - n_true)/max(n_disc, 1)
      }
      discs[[count]] = res
      count = count + 1
    }
  }
  return(list(discs = discs))
}
