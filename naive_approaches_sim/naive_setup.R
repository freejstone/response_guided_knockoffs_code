#!/usr/bin/env Rscript
# Generate the shared designs and responses used by all workers.

source("../functions/functions.R")
source("../functions/helpers.R")
source("../functions/simulation_functions.R")

n = 2000
p = 500
k = 150
nreps = 100
alpha = 0.1
d = 1
delta = 1

set.seed(42)
Sigma = AR1(p = p, a = 3)
X_raw_ar1 = matrix(rnorm(n * p), n, p) %*% chol(Sigma)
X_ar1 = apply(X_raw_ar1, 2, function(x) x / sqrt(sum(x^2)))

set.seed(43)
X_raw_ind = matrix(rnorm(n * p), n, p)
X_ind = apply(X_raw_ind, 2, function(x) x / sqrt(sum(x^2)))

# Optimise ME once for each unique fixed design.
me_optimizer_seed_ar1 = 51001L
me_optimizer_seed_ind = 51002L

cat("Optimising the ME separation vector for the AR1 design...\n")
set.seed(me_optimizer_seed_ar1)
s0_me_ar1 = create_fixed_maxent(
  X_ar1, s0_only = TRUE, verbose = FALSE
)$s0

cat("Optimising the ME separation vector for the independent design...\n")
set.seed(me_optimizer_seed_ind)
s0_me_ind = create_fixed_maxent(
  X_ind, s0_only = TRUE, verbose = FALSE
)$s0

# Scenarios 1 and 3: 90% positive
set.seed(123)
support_1 = sample(1:p, k)
n_neg_1 = round(0.1 * k)
signs_1 = rep(1, k)
signs_1[sample(k, n_neg_1)] = -1
magnitudes_1 = runif(k, min = delta / 2, max = delta)
beta_1 = rep(0, p)
beta_1[support_1] = magnitudes_1 * signs_1

# Scenarios 2 and 4: 90% negative
set.seed(456)
support_2 = sample(1:p, k)
n_pos_2 = round(0.1 * k)
signs_2 = rep(-1, k)
signs_2[sample(k, n_pos_2)] = 1
magnitudes_2 = runif(k, min = delta / 2, max = delta)
beta_2 = rep(0, p)
beta_2[support_2] = magnitudes_2 * signs_2

# Responses
Y_all_ar1_1 = matrix(0, n, nreps)
for (r in 1:nreps) {
  set.seed(1000 + r)
  Y_all_ar1_1[, r] = X_raw_ar1 %*% beta_1 + rnorm(n)
}
Y_all_ar1_2 = matrix(0, n, nreps)
for (r in 1:nreps) {
  set.seed(2000 + r)
  Y_all_ar1_2[, r] = X_raw_ar1 %*% beta_2 + rnorm(n)
}

Y_all_ind_1 = matrix(0, n, nreps)
for (r in 1:nreps) {
  set.seed(3000 + r)
  Y_all_ind_1[, r] = X_raw_ind %*% beta_1 + rnorm(n)
}
Y_all_ind_2 = matrix(0, n, nreps)
for (r in 1:nreps) {
  set.seed(4000 + r)
  Y_all_ind_2[, r] = X_raw_ind %*% beta_2 + rnorm(n)
}

save(X_raw_ar1, X_ar1, X_raw_ind, X_ind,
     s0_me_ar1, s0_me_ind,
     me_optimizer_seed_ar1, me_optimizer_seed_ind,
     beta_1, beta_2,
     Y_all_ar1_1, Y_all_ar1_2, Y_all_ind_1, Y_all_ind_2,
     n, p, k, nreps, alpha, d, delta,
     file = "naive_shared_data.RData")
cat("Saved naive_shared_data.RData\n")
