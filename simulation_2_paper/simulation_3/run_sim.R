#get inputs
args = commandArgs(trailingOnly = TRUE)
seed_ind = as.integer(args[1])
n = as.numeric(args[2])
k = as.numeric(args[3])
if (is.na(seed_ind)) seed_ind = 1
if (is.na(n)) n = 1005
if (is.na(k)) k = 50

print(paste('k, n, seed:', k, n, seed_ind, sep = ' '))

# Libraries
library(hdf5r)
source('../../functions/functions.R')
source('../../functions/simulation_functions.R')
source('../../functions/helpers.R')

# Output directories
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Seed
seed_nk = seed_from_params(n, k, base_seed = 19122025)
set.seed(seed_ind + seed_nk)

# Problem parameters
p = 500     
delta = 0.15
block_size = 5

# Generate the variables from a multivariate normal distribution
Sigma = block_equi_corr(p = 500, block_size = block_size)
groups = rep(seq_len(p / block_size), each = block_size)
X = matrix(rnorm(n*p),n) %*% chol(Sigma)
alphalist = seq(0.3, 0.01, by = -0.01)

# Generate the response from a linear model
beta = create_sparse_coefficients(p = p, sparsity = k/p, 
                                  coeff_dist = "uniform", 
                                  groups = groups,
                                  coeff_size = delta)
y.sample = function(X) X %*% beta + rnorm(n)
y = y.sample(X)
X = apply(X, 2, function(x) x/sqrt(sum(x^2)))

filename = paste("data/sim_",n,'_', k, '_', as.character(seed_ind), ".h5", sep = '')

save_sim_h5(
  filename,
  y, X, beta,
  meta = list(seed = seed_ind, k = k, n = nrow(X), p = ncol(X), model = "linear")
)

# Get noisy response
y_probe = y + rnorm(n, sd = sd(y))

################ sign-guided KFs ################ 

coefs = lm(y_probe ~ X - 1)$coefficients
pos.s = which(coefs > 0)
names(pos.s) = NULL

# We might need to remove more features if n is not large enough for knockoff construction
if (length(pos.s) + p + 2 > n) {
  pos.s = order(coefs, decreasing = TRUE)[1:(n - p - 2)]
  pos.s = sort(pos.s)
}
pos = pos.s

# Response guided KFs
sg_kfs_entropy = response_guided_kf(y, y_probe, X, pos, beta = beta, type = 'entropy', seed = seed_ind + 1052025)

################ Elastic-net-sign-guided KFs ################ 

alphas = c(0.2, 0.4, 0.6, 0.8)
# Use the same CV folds across all alphas so their CV errors are comparable
foldid = sample(rep(seq_len(10), length.out = nrow(X)))
cvs = lapply(alphas, function(a) {
  cv.glmnet(X, y_probe, alpha = a, foldid = foldid, standardize = FALSE, intercept = FALSE)
})

# pick alpha with smallest CV error (at its lambda.min)
mean_cvm = sapply(cvs, function(f) min(f$cvm))
best_i = which.min(mean_cvm)

best_alpha = alphas[best_i]
best_cvfit = cvs[[best_i]]

# pick lambda that minimizes the CV error (lambda.min)
i_min = which.min(best_cvfit$cvm)
cvm_min = best_cvfit$cvm[i_min]

# candidate beta_dense
idx = which(best_cvfit$cvm <= cvm_min)
beta_dense = as.vector(best_cvfit$glmnet.fit$beta[, idx])
pos.s3 = which(beta_dense > 0)

# We might need to remove more features if n is not large enough for knockoff construction
cutoff = 0
while (length(pos.s3) + p + 2 > n) {
  pos.s3 = which(beta_dense > cutoff)
  cutoff = min(beta_dense[beta_dense > cutoff])
}
pos = pos.s3

# Response guided KFs (fall back to OLS sign estimates if enet pos is empty)
if (length(pos) > 0) {
  enet_sg_kfs_entropy = response_guided_kf(y, y_probe, X, pos, beta = beta, type = 'entropy', seed = seed_ind + 1052025)
} else {
  enet_sg_kfs_entropy = sg_kfs_entropy
}

######## Vanilla sign KFs ######## 

kfs_equi = signed_kf(y, X, beta = beta, type = 'equi', seed = seed_ind + 1052025)
kfs_sdp = signed_kf(y, X, beta = beta, type = 'sdp', seed = seed_ind + 1052025)
kfs_me = signed_kf(y, X, beta = beta, type = 'entropy', alphalist = alphalist, seed = seed_ind + 1052025)
kfs_mvr = signed_kf(y, X, beta = beta, type = 'mvr', alphalist = alphalist, seed = seed_ind + 1052025)

savedir = paste('results/results','_',n,'_', k, '_', as.character(seed_ind),'.csv',sep = "")
write.csv(data.frame(SG_entropy_FDP = sg_kfs_entropy$fdp, SG_entropy_power = sg_kfs_entropy$power,

                     enet_SG_entropy_FDP = enet_sg_kfs_entropy$fdp, enet_SG_entropy_power = enet_sg_kfs_entropy$power,

                     equi_FDP = kfs_equi$fdp, equi_power = kfs_equi$power,
                     sdp_FDP = kfs_sdp$fdp, sdp_power = kfs_sdp$power,
                     me_FDP = kfs_me$fdp, me_power = kfs_me$power,
                     mvr_FDP = kfs_mvr$fdp, mvr_power = kfs_mvr$power,
                     
                     alpha = seq(0.3, 0.01, by = -0.01)
                     
), savedir)
