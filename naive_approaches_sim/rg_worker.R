#!/usr/bin/env Rscript
# Run RG(signed, ME) for one scenario and replication.

library(knockoff)
library(glmnet)

source("../functions/functions.R")
source("../functions/helpers.R")
source("../functions/simulation_functions.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript rg_worker.R <task_id>")
task_id = as.integer(args[1])
stopifnot(task_id >= 1, task_id <= 400)

scenario = ceiling(task_id / 100)
rep_id = ((task_id - 1) %% 100) + 1

cat(sprintf("RG task %d: scenario %d, rep %d\n", task_id, scenario, rep_id))

dir.create("results", showWarnings = FALSE)

load("naive_shared_data.RData")

if (scenario == 1) { X = X_ar1; beta = beta_1; y = Y_all_ar1_1[, rep_id]; base_seed = 1000 }
if (scenario == 2) { X = X_ar1; beta = beta_2; y = Y_all_ar1_2[, rep_id]; base_seed = 2000 }
if (scenario == 3) { X = X_ind; beta = beta_1; y = Y_all_ind_1[, rep_id]; base_seed = 3000 }
if (scenario == 4) { X = X_ind; beta = beta_2; y = Y_all_ind_2[, rep_id]; base_seed = 4000 }

seed = base_seed + rep_id
nn = nrow(X)
pp = ncol(X)

set.seed(seed + 30000)
y_probe = y + rnorm(nn, sd = sd(y))
coefs = lm(y_probe ~ X - 1)$coefficients
pos = which(coefs * d > 0)
names(pos) = NULL
if (length(pos) + pp > nn) {
  pos = order(coefs * d, decreasing = TRUE)[1:min(nn - pp, pp)]
  pos = sort(pos)
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

if (length(pos) > 0) {
  res_rg = response_guided_kf(y, y_probe, X, pos, beta = beta, d = d,
                                type = 'entropy', alphalist = c(alpha),
                                seed = seed + 20000)
  fdp_rg = res_rg$fdp[1]
  pow_rg = res_rg$power[1]
} else {
  fdp_rg = 0
  pow_rg = 0
}

out = data.frame(
  scenario = scenario,
  rep = rep_id,
  fdp_rg = fdp_rg,
  pow_rg = pow_rg
)
write.csv(out, sprintf("results/rg_task_%d.csv", task_id), row.names = FALSE)
cat(sprintf("RG task %d done. dFDR: %.3f, Power: %.3f\n",
            task_id, fdp_rg, pow_rg))
