#!/usr/bin/env Rscript
# Run methods (a)-(c) for one scenario and replication.

library(knockoff)
library(glmnet)

source("../functions/functions.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript naive_worker.R <task_id>")
task_id = as.integer(args[1])
stopifnot(task_id >= 1, task_id <= 400)

scenario = ceiling(task_id / 100)
rep_id = ((task_id - 1) %% 100) + 1

cat(sprintf("Naive task %d: scenario %d, rep %d\n", task_id, scenario, rep_id))

dir.create("results", showWarnings = FALSE)

load("naive_shared_data.RData")

if (scenario == 1) { X = X_ar1; beta = beta_1; y = Y_all_ar1_1[, rep_id]; base_seed = 1000 }
if (scenario == 2) { X = X_ar1; beta = beta_2; y = Y_all_ar1_2[, rep_id]; base_seed = 2000 }
if (scenario == 3) { X = X_ind; beta = beta_1; y = Y_all_ind_1[, rep_id]; base_seed = 3000 }
if (scenario == 4) { X = X_ind; beta = beta_2; y = Y_all_ind_2[, rep_id]; base_seed = 4000 }

seed = base_seed + rep_id

compute_W_and_signs = function(X, X_ko, y) {
  p = ncol(X)
  X_all = cbind(X, X_ko)
  swap = sample(1:ncol(X_all))
  Zs = knockoff:::lasso_max_lambda(X_all[, swap], y, intercept = FALSE,
                                     standardize = FALSE, nlambda = 2000)
  Zs[swap] = Zs
  Z = Zs[1:p]
  Zk = Zs[(p + 1):(2 * p)]
  W = pmax(Z, Zk) * sign(Z - Zk)
  W[Z == Zk] = 0
  sign_est = sign(t(X - X_ko) %*% y)
  list(W = W, sign_est = as.vector(sign_est))
}

dir_fdp = function(selected_indices, beta, d) {
  if (length(selected_indices) == 0) return(0)
  n_false = sum(sign(beta[selected_indices]) != d)
  n_false / length(selected_indices)
}

dir_power = function(selected_indices, beta, d) {
  true_pos_d = which(sign(beta) == d)
  if (length(true_pos_d) == 0) return(0)
  sum(selected_indices %in% true_pos_d) / length(true_pos_d)
}

# Methods (a)-(c)
set.seed(seed + 40000)
nn = nrow(X)

ko = knockoff::create.fixed(X, method = "equi", y = y, randomize = TRUE)
X_used = ko$X
y_used = ko$y
X_ko = ko$Xk

res = compute_W_and_signs(X_used, X_ko, y_used)
W = res$W
sign_est = res$sign_est

# (a) Unmodified KF
qvals = TDC_flex_c(W < 0, W > 0, score = abs(W), c = 1/2, lambda = 1/2)
disc_a = which(qvals <= alpha & W > 0)
fdp_a = dir_fdp(disc_a, beta, d)
pow_a = dir_power(disc_a, beta, d)

# (b) Post-hoc sign screening
disc_b_raw = which(qvals <= alpha & W > 0)
disc_b = disc_b_raw[sign_est[disc_b_raw] == d]
fdp_b = dir_fdp(disc_b, beta, d)
pow_b = dir_power(disc_b, beta, d)

# (c) Pre-screening by estimated sign
pre_screen = which(sign_est == d)
if (length(pre_screen) > 0) {
  W_pre = W[pre_screen]
  qvals_pre = TDC_flex_c(W_pre < 0, W_pre > 0, score = abs(W_pre), c = 1/2, lambda = 1/2)
  disc_c_local = which(qvals_pre <= alpha & W_pre > 0)
  disc_c = pre_screen[disc_c_local]
} else {
  disc_c = integer(0)
}
fdp_c = dir_fdp(disc_c, beta, d)
pow_c = dir_power(disc_c, beta, d)

out = data.frame(
  scenario = scenario,
  rep = rep_id,
  fdp_a = fdp_a, fdp_b = fdp_b, fdp_c = fdp_c,
  pow_a = pow_a, pow_b = pow_b, pow_c = pow_c
)
write.csv(out, sprintf("results/naive_task_%d.csv", task_id), row.names = FALSE)
cat(sprintf("Naive task %d done. dFDR: a=%.3f b=%.3f c=%.3f\n",
            task_id, fdp_a, fdp_b, fdp_c))
