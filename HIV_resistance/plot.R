suppressPackageStartupMessages(library(tidyverse))

############ Read in R ############
summarise_R_drug = function(drug, expected_reps = 50) {
  pattern = paste0('^results_[0-9]+_', drug, '[.]csv$')
  files = list.files(path = 'results', pattern = pattern, full.names = TRUE)
  if (length(files) != expected_reps) {
    stop(sprintf('Expected %d R result files for %s, found %d', expected_reps, drug, length(files)))
  }

  required = c('SG_entropy_power', 'enet_SG_entropy_power', 'ME_R_power', 'MVR_R_power',
               'equi_power', 'sdp_power', 'alpha')
  dfs = lapply(files, read.csv)
  missing = vapply(dfs, function(df) any(!required %in% names(df)), logical(1))
  if (any(missing)) {
    stop(sprintf(
      '%d %s result files predate are missing',
      sum(missing), drug
    ))
  }

  bind_rows(dfs) %>%
    group_by(alpha) %>%
    summarise(
      SG_entropy_power = mean(SG_entropy_power),
      enet_SG_entropy_power = mean(enet_SG_entropy_power),
      ME_R_power = mean(ME_R_power),
      MVR_R_power = mean(MVR_R_power),
      equi_power = mean(equi_power),
      sdp_power = mean(sdp_power),
      .groups = 'drop'
    ) %>%
    mutate(sim = drug)
}

df_summary = bind_rows(lapply(c('ATV', 'LPV', 'DRV'), summarise_R_drug))

df_long_summary = df_summary %>%
  pivot_longer(
    cols = -c(alpha, sim),
    names_to = "type",
    values_to = "res"
  )

df_long_summary$sim = factor(
  df_long_summary$sim,
  levels = c("ATV", "LPV", "DRV")
)
inds = sapply(df_long_summary$alpha, function(x) any(abs(x - c(0.05, 0.1, 0.2, 0.3)) <= 1e-5))
df_long_summary = df_long_summary[inds,]
df_long_summary$alpha = factor(df_long_summary$alpha)
df_long_summary = df_long_summary[df_long_summary$type %in% c('SG_entropy_power',
                                                              'enet_SG_entropy_power',
                                                              'ME_R_power',
                                                              'MVR_R_power',
                                                              
                                                              'equi_power',
                                                              'sdp_power'),]
df_long_summary$type[df_long_summary$type == 'SG_entropy_power'] = 'RG(signed, ME)'
df_long_summary$type[df_long_summary$type == 'enet_SG_entropy_power'] = 'RG(reg-signed, ME)'
df_long_summary$type[df_long_summary$type == 'ME_R_power'] = 'S(ME)'
df_long_summary$type[df_long_summary$type == 'MVR_R_power'] = 'S(MVR)'
df_long_summary$type[df_long_summary$type == 'equi_power'] = 'S(equi)'
df_long_summary$type[df_long_summary$type == 'sdp_power'] = 'S(SDP)'

df_long_summary$type = factor(df_long_summary$type, levels = c('RG(signed, ME)',
                                                               'RG(reg-signed, ME)',
                                                               'S(equi)',
                                                               'S(SDP)',
                                                               'S(ME)',
                                                               'S(MVR)'))

set.seed(1234)  # For reproducibility
jitter_offsets <- df_long_summary %>% ungroup() %>%
  distinct(type) %>%
  mutate(
    jitter_x = runif(n(), -0.25, 0.25),
    jitter_y = runif(n(), -0.0, 0.0)
  )

# Apply the jitter to the original data frame
df_jittered <- df_long_summary %>%
  left_join(jitter_offsets, by = c("type")) %>%
  mutate(
    alpha_jittered = as.numeric(alpha) + jitter_x,
    mean_power_jittered = res + jitter_y
  )


p1 <- ggplot(df_jittered) +
  geom_point(aes(x = alpha_jittered, y = mean_power_jittered, colour = type, shape = type), size = 3.5) +
  facet_wrap(. ~ sim, nrow = 1, scales = 'free') +
  scale_shape_manual(values = seq_along(unique(df_long_summary$type))) +
  scale_colour_hue() +
  scale_x_continuous(breaks = seq_along(unique(df_long_summary$alpha)), labels = unique(df_long_summary$alpha)) +
  labs(x = expression(alpha), y = 'Power', shape = 'Method', colour = 'Method') + theme(legend.position = "bottom")



##### Export data #####
df_long_summary$sim_paper = 'HIV'
write.csv(df_long_summary, "HIV_res.csv")

dir.create('plots', showWarnings = FALSE)
ggsave('plots/PI_power.pdf', p1, width = 10, height = 4.5)
