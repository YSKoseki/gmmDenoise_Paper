# 03_histogram.R
# Analyzing the simulation results

# Packages required
library(cowplot); packageVersion("cowplot") # 1.1.1
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Paths
path_in <- "./02-sequencing.R/01-Saved_obj/sequencing.obj"
path_out <- "./03-histogram.R"
dir.create(path_out, recursive = TRUE)
path_obj <- paste0(path_out, "/01-Saved_obj")
dir.create(path_obj, recursive = TRUE)

# Read sequence read data
sequencing <- readRDS(path_in) %>%
  mutate(
    repl = factor(repl),
    seq_type = factor(seq_type, levels = c("True", "False")),
    log_reads = log10(reads)
  )

# Preprocess of data for plotting
xupper <- range(sequencing$log_reads)[2] %>% ceiling()
brk <- seq(0, xupper, .2)
ncol <- length(brk) - 1
preprocess <- function(data = sequencing, breaks = brk) {
  # Get distribution data using hist()
  hist_tib <- data %>%
    nest_by(repl, seq_type) %>%
    mutate(
      hist = list(hist(data$log_reads, breaks = brk, plot = FALSE)),
      mids = list(hist$mids),
      counts = list(hist$counts)
    )
  # "Widen" distribution data table
  colnam <- c("repl", "seq_type", paste0("bin", 1:ncol))
  wide_count_tib <- hist_tib %>%
    unnest(repl) %>% ungroup(repl) %>%
    select(repl, seq_type, counts) %>%
    unnest_wider(counts, names_repair = ~ colnam) %>%
    suppressMessages()
  wide_mid_tib <- hist_tib %>%
    unnest(repl) %>% ungroup(repl) %>%
    select(repl, seq_type, mids) %>%
    unnest_wider(mids, names_repair = ~ colnam) %>%
    suppressMessages()
  # Get mean distribution
  est_count_tib <- wide_count_tib %>%
    select(-repl) %>%
    summarize(across(starts_with("bin"), mean)) %>%
    pivot_longer(cols = starts_with("bin"), values_to = "count")
  mid_tib <- wide_mid_tib %>%
    select(-repl) %>%
    summarize(across(starts_with("bin"), first)) %>%
    pivot_longer(cols = starts_with("bin"), values_to = "mids")
  estimates <- bind_cols(mid_tib, count = est_count_tib$count)
  return(estimates)
}
distr_tib <- preprocess()

# Plot histogram with a zoomed y-axis
plot.hist <- function(
    data = distr_tib, breaks = brk){
  h <- data %>%
    ggplot(aes(x = mids, y = count, fill = seq_type))  +
    geom_bar(width = first(diff(breaks)), stat = "identity") +
    scale_fill_manual(
      values = c(alpha("firebrick", .8), alpha("black", .6)),
      labels = c("True", "False")
    ) +
    scale_x_continuous(limit = c(0, xupper), breaks = seq(0, xupper, 1)) +
    coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    labs(x = "Log10(read size)", y = "Frequency") +
    theme(legend.position = c(.9, .75), legend.justification = c(1, 1)) +
    guides(
      fill = guide_legend(title = NULL),
    )
  return(h)
}
theme_set(cowplot::theme_cowplot())
(hist <- plot.hist())
save_plot(paste0(path_out, "/02-Fig_hist_zoom.svg"), hist, base_asp = 0.9)

# Plot histogram with the full y-axis
plot.hist2 <- function(
    data = distr_tib, breaks = brk){
  h <- data %>%
    ggplot(aes(x = mids, y = count, fill = seq_type))  +
    geom_bar(width = first(diff(breaks)), stat = "identity") +
    scale_fill_manual(
      values = c(alpha("firebrick", .8), alpha("black", .6)),
      labels = c("True", "False")
    ) +
    scale_x_continuous(
      limit = c(0, xupper), breaks = seq(0, xupper, 1), expand = c(0, 0)) +
    scale_y_continuous(
      limit = c(0, 2000), breaks = seq(0, 2000, 500), expand = c(0, 0)
    ) +
    labs(x = "Log10(read size)", y = "Frequency") +
    theme(legend.position = c(.9, .75), legend.justification = c(1, 1)) +
    guides(
      fill = guide_legend(title = NULL),
    )
  return(h)
}
theme_set(cowplot::theme_cowplot())
(hist2 <- plot.hist2())
save_plot(paste0(path_out, "/03-Fig_hist_full.svg"), hist2, base_asp = 0.9)

# A grid plot of the two histograms
theme_set(cowplot::theme_cowplot())
(hist_grid <- cowplot::plot_grid(
  hist2 + theme(axis.title.x=element_blank()),
  hist,
  align = "hv", nrow=2, label_x = c(.87), label_y = c(.95, .95),
  labels=c("(a)", "(b)")
))
save_plot(paste0(path_out, "/04-Fig_histo_grid.svg"), hist_grid,
          base_height = 5, base_asp = .9, ncol = 1, nrow = 1)

# Save R objects
saveRDS(distr_tib, paste0(path_obj, "/distr_tib.obj"))
saveRDS(hist_grid, paste0(path_obj, "/hist_grid.obj"))

# Save workspace and session info
save.image(paste0(path_out, "/histogram.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_out, "/histogram.info"))
