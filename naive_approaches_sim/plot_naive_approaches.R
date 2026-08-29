#!/usr/bin/env Rscript
# Plot empirical directional FDR and power.

library(ggplot2)
library(tidyr)
library(dplyr)

load("naive_approaches_results.RData")

method_cols_fdp = c("fdp_a", "fdp_b", "fdp_c")
method_cols_pow = c("pow_a", "pow_b", "pow_c")
method_names = c("Unmodified KF", "Post-hoc screening", "Pre-screening")

if ("fdp_signed" %in% names(df_all)) {
  method_cols_fdp = c(method_cols_fdp, "fdp_signed")
  method_cols_pow = c(method_cols_pow, "pow_signed")
  method_names = c(method_names, "S(ME)")
}
if ("fdp_rg" %in% names(df_all)) {
  method_cols_fdp = c(method_cols_fdp, "fdp_rg")
  method_cols_pow = c(method_cols_pow, "pow_rg")
  method_names = c(method_names, "RG(signed, ME)")
}

scenario_labels_plot = c(
  "AR1\n90% positive",
  "AR1\n90% negative",
  "Independent\n90% positive",
  "Independent\n90% negative"
)

fdp_long = do.call(rbind, lapply(1:length(method_cols_fdp), function(i) {
  data.frame(
    scenario = scenario_labels_plot[df_all$scenario],
    Method = method_names[i],
    value = df_all[[ method_cols_fdp[i] ]]
  )
}))

pow_long = do.call(rbind, lapply(1:length(method_cols_pow), function(i) {
  data.frame(
    scenario = scenario_labels_plot[df_all$scenario],
    Method = method_names[i],
    value = df_all[[ method_cols_pow[i] ]]
  )
}))

summarise_ci = function(df) {
  df %>%
    group_by(scenario, Method) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      ci_lo = mean_val - qnorm(0.975) * se,
      ci_hi = mean_val + qnorm(0.975) * se,
      .groups = "drop"
    )
}

fdp_summary = summarise_ci(fdp_long)
pow_summary = summarise_ci(pow_long)

scen_levels = scenario_labels_plot
meth_levels = method_names

fdp_summary$scenario = factor(fdp_summary$scenario, levels = scen_levels)
fdp_summary$Method = factor(fdp_summary$Method, levels = meth_levels)
pow_summary$scenario = factor(pow_summary$scenario, levels = scen_levels)
pow_summary$Method = factor(pow_summary$Method, levels = meth_levels)

dodge = position_dodge(width = 0.75)

p_fdr = ggplot(fdp_summary, aes(x = scenario, y = mean_val, fill = Method)) +
  geom_bar(stat = "identity", position = dodge, width = 0.65) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = dodge,
                width = 0.25, linewidth = 0.4) +
  geom_hline(yintercept = alpha, linetype = "dashed", color = "black",
             linewidth = 0.6) +
  scale_fill_hue() +
  labs(x = NULL, y = "Empirical directional FDR", fill = "Method") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank()
  ) +
  annotate("text", x = 0.55, y = alpha + 0.012,
           label = paste0("alpha == ", alpha), parse = TRUE,
           size = 3, hjust = 0)

ggsave("naive_approaches_fdr.pdf", p_fdr, width = 8, height = 5)
cat("Saved naive_approaches_fdr.pdf\n")

p_pow = ggplot(pow_summary, aes(x = scenario, y = mean_val, fill = Method)) +
  geom_bar(stat = "identity", position = dodge, width = 0.65) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = dodge,
                width = 0.25, linewidth = 0.4) +
  scale_fill_hue() +
  labs(x = NULL, y = "Empirical directional Power", fill = "Method") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank()
  )

ggsave("naive_approaches_power.pdf", p_pow, width = 8, height = 5)
cat("Saved naive_approaches_power.pdf\n")
