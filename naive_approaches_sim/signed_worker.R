#!/usr/bin/env Rscript
# Run R-native S(ME) for one scenario and replication.

library(knockoff)
library(glmnet)

source("../functions/functions.R")
source("../functions/helpers.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript signed_worker.R <task_id>")
task_id = as.integer(args[1])
stopifnot(task_id >= 1, task_id <= 400)

scenario = ceiling(task_id / 100)
rep_id = ((task_id - 1) %% 100) + 1

cat(sprintf("Signed task %d: scenario %d, rep %d\n", task_id, scenario, rep_id))

dir.create("results", showWarnings = FALSE)

load("naive_shared_data.RData")

if (!exists("s0_me_ar1") || !exists("s0_me_ind")) {
  stop("Cached ME separation vectors are missing; rerun Rscript naive_setup.R")
}

if (scenario == 1) { X = X_ar1; beta = beta_1; y = Y_all_ar1_1[, rep_id]; base_seed = 1000; s0_me = s0_me_ar1 }
if (scenario == 2) { X = X_ar1; beta = beta_2; y = Y_all_ar1_2[, rep_id]; base_seed = 2000; s0_me = s0_me_ar1 }
if (scenario == 3) { X = X_ind; beta = beta_1; y = Y_all_ind_1[, rep_id]; base_seed = 3000; s0_me = s0_me_ind }
if (scenario == 4) { X = X_ind; beta = beta_2; y = Y_all_ind_2[, rep_id]; base_seed = 4000; s0_me = s0_me_ind }

seed = base_seed + 10000 + rep_id
signed_res = signed_kf(
  y = y,
  X = X,
  beta = beta,
  d = d,
  type = "entropy",
  alphalist = c(alpha),
  seed = seed,
  s0 = s0_me
)
fdp_signed = signed_res$fdp[1]
pow_signed = signed_res$power[1]

out = data.frame(
  scenario = scenario,
  rep = rep_id,
  implementation = "R_native_ME",
  fdp_signed = fdp_signed,
  pow_signed = pow_signed
)
write.csv(out, sprintf("results/signed_task_%d.csv", task_id), row.names = FALSE)
cat(sprintf("Signed task %d done. dFDR: %.3f, Power: %.3f\n",
            task_id, fdp_signed, pow_signed))
