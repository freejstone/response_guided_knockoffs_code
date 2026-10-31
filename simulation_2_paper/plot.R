suppressPackageStartupMessages(library(tidyverse))

args = commandArgs(trailingOnly = TRUE)
preview_mode = any(tolower(args) %in% c("--preview", "preview", "--allow-incomplete"))
since_arg = grep("^--since=", args, value = TRUE)
if (length(since_arg) > 1) {
  stop("Supply at most one --since='YYYY-MM-DD HH:MM:SS' argument")
}
since_time = NULL
if (length(since_arg) == 1) {
  since_text = sub("^--since=", "", since_arg)
  since_time = as.POSIXct(since_text, tz = Sys.timezone())
  if (is.na(since_time)) {
    stop("Could not parse --since; use --since='YYYY-MM-DD HH:MM:SS'")
  }
  preview_mode = TRUE
}
message(if (preview_mode) {
  "Simulation 2 aggregation mode: PREVIEW (incomplete files will be skipped)"
} else {
  "Simulation 2 aggregation mode: STRICT"
})

folders = paste0("simulation_", seq_len(6))
simulation_labels = c(
  simulation_1 = "AR1",
  simulation_2 = "AR1 (Corr)",
  simulation_3 = "Block Equi",
  simulation_4 = "ER (Cov)",
  simulation_5 = "ER (Prec)",
  simulation_6 = "Equi"
)

ns = c(600, 700, 800, 900)
ks = c(50, 100, 150)
seeds = 0:99
alpha_list = seq(0.3, 0.01, by = -0.01)
target_alphas = c(0.05, 0.1, 0.2, 0.3)
# Fixed setting used for the empirical FDR plot in the manuscript.
fdr_n = 700
fdr_k = 100

required_columns = c(
  "SG_entropy_FDP", "SG_entropy_power",
  "enet_SG_entropy_FDP", "enet_SG_entropy_power",
  "equi_FDP", "equi_power",
  "sdp_FDP", "sdp_power",
  "me_FDP", "me_power",
  "mvr_FDP", "mvr_power",
  "alpha"
)

manifest = crossing(folder = folders, n = ns, k = ks, seed = seeds) %>%
  mutate(path = file.path(folder, "results", sprintf("results_%d_%d_%d.csv", n, k, seed)))

missing_files = manifest$path[!file.exists(manifest$path)]
if (length(missing_files) > 0) {
  preview = paste(head(missing_files, 20), collapse = "\n  ")
  problem = sprintf(
    "Missing %d of %d expected R result CSV files. First missing files:\n  %s",
    length(missing_files), nrow(manifest), preview
  )
  if (preview_mode) warning(problem, immediate. = TRUE) else stop(problem)
}

manifest_to_read = manifest %>% filter(file.exists(path))
if (!is.null(since_time)) {
  modification_times = file.info(manifest_to_read$path)$mtime
  manifest_to_read = manifest_to_read[modification_times >= since_time, , drop = FALSE]
}
if (nrow(manifest_to_read) == 0) {
  stop("No eligible result files were found")
}

read_result = function(folder, n, k, seed, path) {
  reject_file = function(reason) {
    problem = sprintf("Skipping %s: %s", path, reason)
    if (preview_mode) {
      warning(problem, immediate. = TRUE)
      return(NULL)
    }
    stop(problem)
  }

  size = file.info(path)$size
  if (is.na(size) || size == 0) {
    return(reject_file("file is empty or unavailable"))
  }

  read_error = NULL
  df = tryCatch(
    read.csv(path, check.names = FALSE),
    error = function(e) {
      read_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(df)) {
    return(reject_file(paste("could not be read:", read_error)))
  }

  missing_columns = setdiff(required_columns, names(df))
  if (length(missing_columns) > 0) {
    return(reject_file(sprintf(
      "missing required active R columns: %s",
      paste(missing_columns, collapse = ", ")
    )))
  }

  if (nrow(df) != length(alpha_list)) {
    return(reject_file(sprintf(
      "has %d rows; expected %d alpha values",
      nrow(df), length(alpha_list)
    )))
  }
  observed_alpha = sort(as.numeric(df$alpha))
  expected_alpha = sort(alpha_list)
  if (any(!is.finite(observed_alpha)) ||
      !isTRUE(all.equal(observed_alpha, expected_alpha, tolerance = 1e-10,
                        check.attributes = FALSE))) {
    return(reject_file("does not contain the expected alpha grid"))
  }

  active_columns = setdiff(required_columns, "alpha")
  active_values = as.matrix(df[, active_columns, drop = FALSE])
  if (!is.numeric(active_values) || any(!is.finite(active_values))) {
    return(reject_file("contains non-finite or non-numeric active results"))
  }

  df %>%
    select(all_of(required_columns)) %>%
    mutate(folder = folder, n = n, k = k, seed = seed)
}

raw_results = pmap_dfr(manifest_to_read, read_result)
if (nrow(raw_results) == 0) {
  stop("None of the eligible result files was complete and readable")
}
if (preview_mode) {
  message(sprintf(
    "Preview: using %d complete result files out of %d expected for the final analysis",
    nrow(raw_results) / length(alpha_list), nrow(manifest)
  ))
}

power_labels = c(
  SG_entropy_power = "RG(signed, ME)",
  enet_SG_entropy_power = "RG(reg-signed, ME)",
  equi_power = "S(equi)",
  sdp_power = "S(SDP)",
  me_power = "S(ME)",
  mvr_power = "S(MVR)"
)
fdp_labels = c(
  SG_entropy_FDP = "RG(signed, ME)",
  enet_SG_entropy_FDP = "RG(reg-signed, ME)",
  equi_FDP = "S(equi)",
  sdp_FDP = "S(SDP)",
  me_FDP = "S(ME)",
  mvr_FDP = "S(MVR)"
)
method_levels = unique(unname(power_labels))

summarise_results_by_setting = function(df, labels) {
  df %>%
    pivot_longer(
      cols = all_of(names(labels)),
      names_to = "source_type",
      values_to = "res"
    ) %>%
    mutate(
      alpha = vapply(alpha, function(value) {
        distance = abs(target_alphas - value)
        if (min(distance) <= 1e-8) target_alphas[which.min(distance)] else NA_real_
      }, numeric(1)),
      type = unname(labels[source_type])
    ) %>%
    filter(!is.na(alpha)) %>%
    # Average only over Monte Carlo repetitions. Retain n and k so that
    # relative power can later be calculated within each exact setting.
    group_by(folder, n, k, alpha, type) %>%
    summarise(res = mean(res), .groups = "drop") %>%
    mutate(
      sim = factor(unname(simulation_labels[folder]), levels = unname(simulation_labels)),
      type = factor(type, levels = method_levels)
    ) %>%
    select(sim, n, k, alpha, type, res)
}

average_over_design_grid = function(df) {
  df %>%
    group_by(sim, alpha, type) %>%
    summarise(res = mean(res), .groups = "drop")
}

power_summary_by_setting = summarise_results_by_setting(raw_results, power_labels)
fdp_summary_by_setting = summarise_results_by_setting(raw_results, fdp_labels)

# Power is averaged over the (n, k) grid; FDR uses only the manuscript setting.
power_summary = average_over_design_grid(power_summary_by_setting)
fdp_summary = fdp_summary_by_setting %>% filter(n == fdr_n, k == fdr_k)
if (nrow(fdp_summary) == 0) {
  stop(sprintf("No eligible results for the empirical FDR setting n = %d, k = %d", fdr_n, fdr_k))
}

make_jittered = function(df) {
  offsets = df %>%
    distinct(type) %>%
    mutate(jitter_x = runif(n(), -0.25, 0.25))

  df %>%
    left_join(offsets, by = "type") %>%
    mutate(alpha_x = match(alpha, target_alphas) + jitter_x)
}

set.seed(1234)
power_jittered = make_jittered(power_summary)
shape_values = setNames(seq_along(method_levels), method_levels)

p1 = ggplot(power_jittered) +
  geom_point(aes(x = alpha_x, y = res, colour = type, shape = type), size = 3.5) +
  facet_wrap(~ sim, nrow = 2, scales = "free_y") +
  scale_shape_manual(values = shape_values) +
  scale_colour_hue() +
  scale_x_continuous(breaks = seq_along(target_alphas), labels = target_alphas) +
  labs(x = expression(alpha), y = "Power", shape = "Method", colour = "Method") +
  theme(legend.position = "bottom")

set.seed(1234)
fdp_jittered = make_jittered(fdp_summary)
ref_df = data.frame(
  x = seq_along(target_alphas) - 0.4,
  xend = seq_along(target_alphas) + 0.4,
  y = target_alphas
)

p2 = ggplot(fdp_jittered) +
  geom_segment(
    data = ref_df,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_point(aes(x = alpha_x, y = res, colour = type, shape = type), size = 3.5) +
  facet_wrap(~ sim, nrow = 2, scales = "free_y") +
  scale_shape_manual(values = shape_values) +
  scale_colour_hue() +
  scale_x_continuous(breaks = seq_along(target_alphas), labels = target_alphas) +
  labs(x = expression(alpha), y = "FDR", shape = "Method", colour = "Method") +
  theme(legend.position = "bottom")

# Export the setting-level powers for within-(model, n, k, alpha)
# normalisation in code/plot_ratio.R.
df_long_summary = power_summary_by_setting %>% mutate(sim_paper = "simulation_2")
write.csv(df_long_summary, "sim_2_res.csv", row.names = FALSE)

dir.create("plots", recursive = TRUE, showWarnings = FALSE)
ggsave("plots/simulation_2_power.pdf", p1, width = 10, height = 7)
ggsave("plots/simulation_2_fdr.pdf", p2, width = 10, height = 7)
