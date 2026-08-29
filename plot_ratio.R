suppressMessages(library(tidyverse))

base = "/Users/jfre0619/2026_response_guided_knockoffs/code"

files = c("Simulation 1"   = "simulation_1_paper/sim_1_res.csv",
           "Simulation 2"   = "simulation_2_paper/sim_2_res.csv",
           "HIV resistance" = "HIV_resistance/HIV_res.csv")

method_levels = c("RG(signed, ME)", "RG(reg-signed, ME)",
                   "S(equi)", "S(SDP)", "S(ME)", "S(MVR)")

# default ggplot hue palette used by the existing plot.R figures (scale_colour_hue),
method_cols = c("RG(signed, ME)"     = "#F8766D",
                 "RG(reg-signed, ME)" = "#B79F00",
                 "S(equi)"            = "#00BA38",
                 "S(SDP)"             = "#00BFC4",
                 "S(ME)"              = "#619CFF",
                 "S(MVR)"             = "#F564E3")

dat = imap_dfr(files, function(path, lab) {
  df = read.csv(file.path(base, path), stringsAsFactors = FALSE)

  # The simulations retain their n and k settings. The HIV analysis has no
  # analogous design grid, so NA placeholders keep each (drug, alpha) cell.
  if (lab == "HIV resistance") {
    df$n = NA_real_
    df$k = NA_real_
  } else if (!all(c("n", "k") %in% names(df))) {
    stop(sprintf(
      "%s does not contain n and k; rerun that simulation's plot.R first",
      file.path(base, path)
    ))
  }

  df %>%
    transmute(dataset = lab, group = sim, n, k, alpha, type, res)
}) %>%
  filter(type %in% method_levels) %>%
  mutate(dataset = factor(dataset, levels = names(files)),
         type    = factor(type, levels = method_levels))

setting_counts = dat %>%
  distinct(dataset, group, n, k, alpha) %>%
  count(dataset, name = "settings")
expected_settings = c(
  "Simulation 1" = 360L,
  "Simulation 2" = 288L,
  "HIV resistance" = 12L
)
for (i in seq_len(nrow(setting_counts))) {
  dataset_name = as.character(setting_counts$dataset[i])
  observed = setting_counts$settings[i]
  expected = unname(expected_settings[dataset_name])
  message(sprintf(
    "%s: using %d of %d expected relative-power settings",
    dataset_name, observed, expected
  ))
  if (observed < expected) {
    warning(sprintf(
      "%s is incomplete; this figure is a preview and must not be used as the final result",
      dataset_name
    ), immediate. = TRUE)
  }
}

# Compare methods only within the same data-generating model, n, k, and alpha.
# For HIV, n and k are NA and the comparison remains within (drug, alpha).
ratio_df = dat %>%
  group_by(dataset, group, n, k, alpha) %>%
  mutate(
    best_power = max(res),
    # If every method has zero power, all methods are tied for best.
    ratio = if_else(best_power > 0, res / best_power, 1)
  ) %>%
  ungroup()

row1 = ratio_df %>%
  transmute(dataset, type, ratio, panel_row = "Relative Power")

# Row 2 averages the setting-level ratios shown in row 1: 360 for Simulation 1,
# 288 for Simulation 2, and 12 for HIV, per method.
row2 = ratio_df %>%
  group_by(dataset, type) %>%
  summarise(ratio = mean(ratio), .groups = "drop") %>%
  mutate(panel_row = "Mean Relative Power")

panel_levels = c("Relative Power", "Mean Relative Power")
row1$panel_row = factor(row1$panel_row, levels = panel_levels)
row2$panel_row = factor(row2$panel_row, levels = panel_levels)

make_ratio_plot = function(relative_data, mean_data) {
  ggplot(mapping = aes(x = type, y = ratio, fill = type)) +
    geom_boxplot(data = relative_data, outlier.shape = NA, alpha = 0.55,
                 width = 0.65, linewidth = 0.3) +
    geom_jitter(data = relative_data, width = 0.12, height = 0,
                size = 0.7, alpha = 0.5) +
    geom_col(data = mean_data, alpha = 0.85, width = 0.7) +
    facet_grid(panel_row ~ dataset) +
    scale_fill_manual(values = method_cols, guide = "none") +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          panel.grid.minor = element_blank())
}

make_standalone_ratio_plot = function(relative_data, mean_data, dataset_name) {
  ggplot(mapping = aes(x = type, y = ratio, fill = type)) +
    geom_boxplot(data = relative_data, outlier.shape = NA, alpha = 0.55,
                 width = 0.65, linewidth = 0.3) +
    geom_jitter(data = relative_data, width = 0.12, height = 0,
                size = 0.7, alpha = 0.5) +
    geom_col(data = mean_data, alpha = 0.85, width = 0.7) +
    facet_grid(. ~ panel_row) +
    scale_fill_manual(values = method_cols, guide = "none") +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
    labs(title = dataset_name, x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

p = make_ratio_plot(row1, row2)

dir.create(file.path(base, "summary_plots"), showWarnings = FALSE)
ggsave(file.path(base, "summary_plots", "tpr_ratio_2x3.pdf"), p,
       width = 11, height = 7, device = cairo_pdf)
ggsave(file.path(base, "summary_plots", "tpr_ratio_2x3.png"), p,
       width = 11, height = 7, dpi = 150)

standalone_files = c(
  "Simulation 1" = "tpr_ratio_simulation_1.pdf",
  "Simulation 2" = "tpr_ratio_simulation_2.pdf",
  "HIV resistance" = "tpr_ratio_hiv.pdf"
)
for (dataset_name in names(standalone_files)) {
  standalone_plot = make_standalone_ratio_plot(
    row1 %>% filter(dataset == dataset_name),
    row2 %>% filter(dataset == dataset_name),
    dataset_name
  )
  ggsave(
    file.path(base, "summary_plots", standalone_files[[dataset_name]]),
    standalone_plot,
    width = 12,
    height = 5,
    device = cairo_pdf
  )
}

message(paste(
  "Saved combined summary_plots/tpr_ratio_2x3.{pdf,png}",
  "and three landscape standalone PDF panels"
))
