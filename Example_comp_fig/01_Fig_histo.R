# 01_Fig_histo.R
# Last updated on 2023.1.4 by YK
# An R script to plot a composite histogram plot
# R 4.1.2

# Packages required
library(cowplot); packageVersion("cowplot") # 1.1.1
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Load objects
fig_histo_st <- readRDS("../Example_1/11-GMMdn_DRA005106/01-Saved_object/fig_histo.obj")
fig_histo_es <- readRDS("../Example_2/11-GMMdn_DRA010703/01-Saved_object/fig_histo.obj")

# Output directory
path_output <- "./01-Fig_histo"
dir.create(path_output, recursive=TRUE)

#The composite histogram plot
theme_set(cowplot::theme_cowplot())
(fig_histo_grid <- cowplot::plot_grid(fig_histo_st[["uno2"]] +
                                        scale_x_continuous(limits = c(1,6)) +
                                        scale_y_continuous(
                                          expand = c(0,0), limits = c(0,150)
                                        ),
                                      fig_histo_st[["uno5"]] +
                                        scale_x_continuous(limits = c(1,6)) +
                                        scale_y_continuous(
                                          expand = c(0,0), limits = c(0,150)
                                        ),
                                      fig_histo_es[["uno2"]] +
                                        scale_x_continuous(limits = c(1,6)) +
                                        scale_y_continuous(
                                          expand = c(0,0), limits = c(0,150)
                                        ),
                                      fig_histo_es[["uno5"]] +
                                        scale_x_continuous(limits = c(1,6)) +
                                        scale_y_continuous(
                                          expand = c(0,0), limits = c(0,150)
                                        ),
                                      labels=c(
                                        "(a) Stream, \u03b1 = 2",
                                        "(b) Stream, \u03b1 = 5",
                                        "(c) Estuarine, \u03b1 = 2",
                                        "(d) Estuarine, \u03b1 = 5"
                                      ),
                                      label_x=.18, label_y=1, hjust=-.1, vjust=2.2)
)
save_plot(paste0(path_output, "/01-Fig_histo_grid.svg"),
          fig_histo_grid,
          ncol=2,
          nrow=2,
          base_asp = 1
)

# Save workspace and session info
save.image(paste0(path_output, "/fig_histo.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/fig_histo.info"))
