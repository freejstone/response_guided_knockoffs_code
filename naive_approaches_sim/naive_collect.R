#!/usr/bin/env Rscript
# Usage: Rscript naive_collect.R [naive|all]

args = commandArgs(trailingOnly = TRUE)
mode = if (length(args) >= 1) args[1] else "all"

alpha = 0.1
scenario_labels = c(
  "AR1, 90% positive",
  "AR1, 90% negative",
  "Independent, 90% positive",
  "Independent, 90% negative"
)

validate_task_grid = function(df, label) {
  expected = expand.grid(scenario = 1:4, rep = 1:100)
  expected_keys = paste(expected$scenario, expected$rep, sep = ":")
  observed_keys = paste(df$scenario, df$rep, sep = ":")
  missing_keys = setdiff(expected_keys, observed_keys)
  duplicated_keys = unique(observed_keys[duplicated(observed_keys)])
  unexpected_keys = setdiff(observed_keys, expected_keys)

  if (length(missing_keys) > 0 || length(duplicated_keys) > 0 ||
      length(unexpected_keys) > 0) {
    stop(sprintf(
      "%s results are incomplete or duplicated: %d missing, %d duplicated, %d unexpected task keys",
      label, length(missing_keys), length(duplicated_keys), length(unexpected_keys)
    ))
  }
  invisible(TRUE)
}

# Naive results
naive_files = list.files("results", pattern = "^naive_task_[0-9]+\\.csv$",
                          full.names = TRUE)
cat(sprintf("Found %d naive result files\n", length(naive_files)))

if (length(naive_files) == 0) stop("No naive result files found in results/")
df_naive = do.call(rbind, lapply(naive_files, read.csv))
validate_task_grid(df_naive, "Naive")

if (mode == "naive") {
  for (s in 1:4) {
    sub = df_naive[df_naive$scenario == s, ]
    if (nrow(sub) == 0) next
    cat(sprintf("\n%s (%d reps) — target alpha = %g:\n",
                scenario_labels[s], nrow(sub), alpha))
    cat(sprintf("  (a) Unmodified KF:           dFDR=%.3f  Power=%.3f\n",
                mean(sub$fdp_a), mean(sub$pow_a)))
    cat(sprintf("  (b) Post-hoc sign screening: dFDR=%.3f  Power=%.3f\n",
                mean(sub$fdp_b), mean(sub$pow_b)))
    cat(sprintf("  (c) Pre-screening by sign:   dFDR=%.3f  Power=%.3f\n",
                mean(sub$fdp_c), mean(sub$pow_c)))
  }
  cat("\n")
  quit("no")
}

# Signed results
signed_files = list.files("results", pattern = "^signed_task_[0-9]+\\.csv$",
                           full.names = TRUE)
cat(sprintf("Found %d signed result files\n", length(signed_files)))
df_signed = if (length(signed_files) > 0) {
  signed_rows = lapply(signed_files, read.csv)
  valid_implementation = vapply(
    signed_rows,
    function(x) {
      "implementation" %in% names(x) &&
        nrow(x) == 1 &&
        identical(as.character(x$implementation), "R_native_ME")
    },
    logical(1)
  )
  if (!all(valid_implementation)) {
    stop(paste(
      "Signed results contain legacy or unknown implementations.",
      "Rerun every signed task with the R-native signed_worker.R before aggregation."
    ))
  }
  ans = do.call(rbind, signed_rows)
  validate_task_grid(ans, "Signed")
  ans
} else { NULL }
if (mode == "all" && is.null(df_signed)) {
  stop("No R-native signed S(ME) result files found in results/")
}

# RG results
rg_files = list.files("results", pattern = "^rg_task_[0-9]+\\.csv$",
                       full.names = TRUE)
cat(sprintf("Found %d RG result files\n", length(rg_files)))
df_rg = if (length(rg_files) > 0) {
  ans = do.call(rbind, lapply(rg_files, read.csv))
  validate_task_grid(ans, "RG")
  ans
} else { NULL }
if (mode == "all" && is.null(df_rg)) {
  stop("No response-guided result files found in results/")
}

df_all = df_naive
if (!is.null(df_signed)) {
  df_all = merge(df_all, df_signed, by = c("scenario", "rep"), all.x = TRUE)
}
if (!is.null(df_rg)) {
  df_all = merge(df_all, df_rg, by = c("scenario", "rep"), all.x = TRUE)
}

for (s in 1:4) {
  sub = df_all[df_all$scenario == s, ]
  if (nrow(sub) == 0) next
  cat(sprintf("\n%s (%d reps) — target alpha = %g:\n",
              scenario_labels[s], nrow(sub), alpha))
  cat(sprintf("  (a) Unmodified KF:           dFDR=%.3f  Power=%.3f\n",
              mean(sub$fdp_a), mean(sub$pow_a)))
  cat(sprintf("  (b) Post-hoc sign screening: dFDR=%.3f  Power=%.3f\n",
              mean(sub$fdp_b), mean(sub$pow_b)))
  cat(sprintf("  (c) Pre-screening by sign:   dFDR=%.3f  Power=%.3f\n",
              mean(sub$fdp_c), mean(sub$pow_c)))
  if ("fdp_signed" %in% names(sub))
    cat(sprintf("  (d) Signed KF S(ME):         dFDR=%.3f  Power=%.3f\n",
                mean(sub$fdp_signed, na.rm = TRUE),
                mean(sub$pow_signed, na.rm = TRUE)))
  if ("fdp_rg" %in% names(sub))
    cat(sprintf("  (e) RG(signed, ME):          dFDR=%.3f  Power=%.3f\n",
                mean(sub$fdp_rg, na.rm = TRUE),
                mean(sub$pow_rg, na.rm = TRUE)))
}

load("naive_shared_data.RData")
save(df_all, beta_1, beta_2, alpha, scenario_labels,
     file = "naive_approaches_results.RData")
cat("\nSaved naive_approaches_results.RData\n")
