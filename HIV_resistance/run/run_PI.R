###### libraries ######
library(tidyverse)
source('../../functions/functions.R')
source('../../functions/helpers.R')

###### Parameters ###### 

alphalist = seq(0.3, 0.01, by = -0.01)
reps = 50

###### Seed ######
seed = 6122025
set.seed(seed)

###### Output directory ######
dir.create("../results", showWarnings = FALSE)

###### Process ###### 

drug_class = 'PI'

gene_df = read.delim('RB_data/PI_DataSet.txt', na.string = c('NA', ''), stringsAsFactors = FALSE)
truth_df = read.csv('RB_data/PI_truth.csv')
position_mutations = truth_df$Condition[grepl("^\\d+[A-Z]+$", truth_df$Condition)]

# Returns rows for which every column matches the given regular expression.
grepl_rows <- function(pattern, df) {
  cell_matches = apply(df, c(1,2), function(x) grepl(pattern, x))
  apply(cell_matches, 1, all)
}

pos_start = which(names(gene_df) == 'P1')
pos_cols = seq.int(pos_start, ncol(gene_df) - 1)
valid_rows = grepl_rows('^(\\.|-|[A-Zid]+)$', gene_df[,pos_cols])
gene_df = gene_df[valid_rows,]

# Flatten a matrix to a vector with names from concatenating row/column names.
flatten_matrix <- function(M, sep='.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep=sep)))
  x
}

# Construct preliminary design matrix.
muts = c(LETTERS, 'i', 'd')
X = outer(muts, as.matrix(gene_df[,pos_cols]), Vectorize(grepl))
X = aperm(X, c(2,3,1))
dimnames(X)[[3]] <- muts
X = t(apply(X, 1, flatten_matrix))
mode(X) <- 'numeric'

# Remove any mutation/position pairs that never appear in the data.
X = X[,colSums(X) != 0]

# Extract response matrix.
Y = gene_df[, c("ATV", "LPV", "DRV")] 

process <- function (X, y) {
  # Log-transform the drug resistance measurements.
  y = log(y)
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  X = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  X = X[,colSums(X) >= 3]
  
  # Remove duplicate predictors.
  X = X[,colSums(abs(cor(X)-1) < 1e-4) == 1]

  # Use the unit-norm convention assumed by the fixed-X methods and paper.
  X = sweep(X, 2, sqrt(colSums(X^2)), "/")
  
  list(X = X, y = y)
}

results = lapply(Y, function(y) process(X, y))

negative_signs = data.frame(
  Position = c(50, 50, 82, 82, 82, 88),
  Mutation = c("L", "L", "A", "T", "S", "S"),
  type = c("LPV", "DRV", "DRV", "DRV", "DRV", "DRV")
)

################ Doing KFs ################ 

for (j in 1:length(results)) {
  drug.name = names(results)[j]
  y = results[[j]]$y
  X = results[[j]]$X
  
  negative.drug = negative_signs[negative_signs$type == drug.name, ]
  
  ## extract numeric position
  pos <- as.integer(sub("(\\d+).*", "\\1", position_mutations))
  ## extract letter part
  mut_str <- sub("\\d+", "", position_mutations)
  ## split letters and build data frame
  df.pos_mut <- data.frame(
    Position = rep(pos, nchar(mut_str)),
    Mutation = unlist(strsplit(mut_str, "")),
    row.names = NULL
  )
  
  df.pos_mut <- anti_join(df.pos_mut, negative.drug, by = c("Position", "Mutation"))
  pos_mut <- paste0(df.pos_mut$Position, df.pos_mut$Mutation)
  
  # Column names look like "P74.P" (position 74, mutation letter P). Strip only the
  # leading "P" and the "." so the mutation letter is kept, giving "74P".
  stopifnot(all(grepl("^P\\d+\\.[A-Za-z]$", colnames(X))))
  position_mutation_X = sub("^P(\\d+)\\.(.)$", "\\1\\2", colnames(X))
  
  beta = rep(0, ncol(X))
  beta[position_mutation_X %in% pos_mut] = 1
  
  n = nrow(X)
  p = ncol(X)

  # All four full-design separation vectors depend only on X, so compute them once per drug.
  # with_seed restores the existing RNG stream used for the response probes.
  ME_R_s0 = withr::with_seed(
    seed + 100000 + j,
    create_fixed_maxent(X, s0_only = TRUE, verbose = FALSE)$s0
  )
  MVR_R_s0 = withr::with_seed(
    seed + 200000 + j,
    create_fixed_mvr(X, s0_only = TRUE, verbose = FALSE)$s0
  )
  equi_s0 = withr::with_seed(
    seed + 300000 + j,
    create_fixed_equi_sdp(X, type = 'equi', s0_only = TRUE, verbose = FALSE)$s0
  )
  sdp_s0 = withr::with_seed(
    seed + 400000 + j,
    create_fixed_equi_sdp(X, type = 'sdp', s0_only = TRUE, verbose = FALSE)$s0
  )
  
  for (rep in 1:reps) {
    # Start the probing
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
    sg_kfs_entropy = response_guided_kf(y, y_probe, X, pos, beta = beta, type = 'entropy', seed = rep + seed)
    
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
    
    # Response guided KFs
    enet_sg_kfs_entropy = response_guided_kf(y, y_probe, X, pos, beta = beta, type = 'entropy', seed = rep + seed)
    
    ######## Standard signed KFs ########
    # Reuse each optimized separation vector, while redrawing the randomized
    # knockoff complement and refitting the LSM scorer in every replicate.
    kfs_me_R = signed_kf(y, X, beta = beta, type = 'entropy', seed = rep + seed, s0 = ME_R_s0)
    kfs_mvr_R = signed_kf(y, X, beta = beta, type = 'mvr', seed = rep + seed, s0 = MVR_R_s0)
    kfs_equi = signed_kf(y, X, beta = beta, type = 'equi', seed = rep + seed, s0 = equi_s0)
    kfs_sdp = signed_kf(y, X, beta = beta, type = 'sdp', seed = rep + seed, s0 = sdp_s0)
    
    savedir = paste('../results/results','_',rep,'_', drug.name, '.csv',sep = "")
    write.csv(data.frame(
      SG_entropy_FDP = sg_kfs_entropy$fdp, SG_entropy_power = sg_kfs_entropy$power*sum(beta),

      enet_SG_entropy_FDP = enet_sg_kfs_entropy$fdp, enet_SG_entropy_power = enet_sg_kfs_entropy$power*sum(beta),

      ME_R_FDP = kfs_me_R$fdp, ME_R_power = kfs_me_R$power*sum(beta),
      MVR_R_FDP = kfs_mvr_R$fdp, MVR_R_power = kfs_mvr_R$power*sum(beta),
      equi_FDP = kfs_equi$fdp, equi_power = kfs_equi$power*sum(beta),
      sdp_FDP = kfs_sdp$fdp, sdp_power = kfs_sdp$power*sum(beta),
      
      alpha = seq(0.3, 0.01, by = -0.01)
      
    ), savedir)
  }
}
