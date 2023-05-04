# 11_GMM_DRA010703.R
# Last updated on 2023.5.1 by YK
# An R script to infer read count cut-off threshold for ASVs, based on Gaussian mixture modeling (GMM)
# R 4.1.2

# Packages required
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9018
library(mixtools); packageVersion("mixtools") # 1.2.0
library(cowplot); packageVersion("cowplot") # 1.1.1
library(wesanderson); packageVersion("wesanderson") # 0.3.6
library(tidyverse); packageVersion("tidyverse") # 1.3.1
library(gmmDenoise); packageVersion("gmmDenoise") # 0.2.3

# Paths
## Input: the list of phyloseq objects
path_input <- "./10-Phyloseq_DRA010703/01-Saved_object/phylo.obj"
## Output directory
path_output <- "./11-GMMdn_DRA010703"

# Load the phyloseq objects
phylo <- readRDS(path_input)
reads <- phylo %>%
  lapply(function(x) {x %>% otu_table() %>% colSums()})

# Step 0: Plot ASV read count histogram (in log scale)
fig_histo <- reads %>% 
  lapply(function(x) asvhist(x, scale="log", xlim = c(1, 6), ylim=c(0, 500)))
theme_set(cowplot::theme_cowplot())
(fig_histo_grid <- fig_histo[c(4, 1:3)] %>% 
  cowplot::plot_grid(plotlist=., 
                     labels=toupper(names(fig_histo)[c(4, 1:3)]), 
                     label_x=1, label_y=1, hjust=2, vjust=2))
dir.create(path_output, recursive=TRUE)
save_plot(paste0(path_output, "/02-Fig_histo_grid.svg"), fig_histo_grid, ncol=2, nrow=2)

# Step 1: Selecting the number of components
## Approach 1: Cross-validation
### The CV analysis
set.seed(101)
crossval <- reads %>% 
  lapply(function(x) gmmcv(log10(x), maxit=5000, epsilon=1e-2))
### Plot the CV results
fig_cv <- crossval %>% lapply(autoplot)
theme_set(cowplot::theme_cowplot())
(fig_cv_grid <- fig_cv[c(4, 1:3)] %>% 
    cowplot::plot_grid(plotlist=., 
                       labels=toupper(names(fig_cv)[c(4, 1:3)]), 
                       label_x=1, label_y=1, hjust=2, vjust=2))
save_plot(paste0(path_output, "/03-Fig_cv_grid.svg"), fig_cv_grid, ncol=2, nrow=2)

# Step 2: Model fitting with the selected number of components
## The selected number of components
k <- list(data=2, uno5=3, uno2=3, mifs=3)
## Fit k-component GMMs
set.seed(103)
gmm_fit <- map2(reads, k, 
                function(x, k) gmmem(log10(x), k, maxit=5000, epsilon=1e-2))

# Step 3: Plot GMM-estimated count density functions (CDF) on histogram
## Draw the plots
fig_pdf <- gmm_fit %>%
  lapply(function(x) {
    autoplot(x, xlim = c(1, 6), ylim=c(0, 500))
  })
theme_set(cowplot::theme_cowplot())
fig_pdf_grid <- fig_pdf[c(4, 1:3)] %>% 
  cowplot::plot_grid(plotlist=., 
                     labels=toupper(names(fig_pdf)[c(4, 1:3)]), 
                     label_x=.75, label_y=1, hjust=0, vjust=2)
save_plot(paste0(path_output, "/04-Fig_pdf_grid.svg"), fig_pdf_grid, ncol=2, nrow=2)

# Step 4: Infer read count threshold for ASV cut-off
## Calculate 95-percentiles of individual mixture components of fitted GMMs
lower95_log <- gmm_fit %>% 
  lapply(function(x) {
    ncomp <- x$mu %>% length()
    quantile(x, comp=1:ncomp)
  })
(thresh_tab <- list(dada=c(1, NA), uno5=c(NA, 1, NA), uno2=c(NA, 1, NA), mifs=c(NA, 1, NA)) %>% 
    map2(lower95_log, function(x, y) x*y) %>% 
    lapply(function(x) {
      if (all(is.na(x))) NA
      else na.omit(x)
    }) %>% 
    unlist() %>% 
    data.frame(log=.) %>% 
    rownames_to_column(var="data") %>% 
    mutate(norm=ceiling(10^log)))
write.csv(thresh_tab[c(4, 1:3), ], paste0(path_output, "/05-thresh_tab.csv"))

## Draw vertical lines indicating 95-percentiles to the plots
fig_pdf2 <- list(dada=c(1, NA), uno5=c(NA, 1, NA), uno2=c(NA, 1, NA), mifs=c(NA, 1, NA)) %>% 
  map2(lower95_log, function(x, y) x*y) %>% 
  map2(gmm_fit, .,
       function(x, y) {
         autoplot(x, vline=y, xlim = c(1, 6), ylim=c(0, 500))
       })
theme_set(cowplot::theme_cowplot())
fig_pdf_grid2 <- fig_pdf2[c(4, 1:3)] %>% 
    cowplot::plot_grid(plotlist=., 
                       labels=toupper(names(fig_pdf)[c(4, 1:3)]), 
                       label_x=.75, label_y=1, hjust=0, vjust=2)
save_plot(paste0(path_output, "/06-Fig_pdf_grid2.svg"), fig_pdf_grid2, ncol=2, nrow=2)

# Compound plot for publication
theme_set(cowplot::theme_cowplot())
(fig_publ <- cowplot::plot_grid(
  fig_histo[["mifs"]] +
    scale_x_continuous(limits=c(1, 6), breaks=seq(1, 6)) +
    scale_y_continuous(limits=c(0, 510), breaks=seq(0, 500, 100), expand=c(0, 0)) +
    xlab(""),
  fig_histo[["dada"]] +
    scale_x_continuous(limits=c(1, 6), breaks=seq(1, 6)) +
    scale_y_continuous(limits=c(0, 155), breaks=seq(0, 150, 50), expand=c(0, 0)) +
    theme(axis.title.y=element_blank()),
  fig_histo[["uno5"]] +
    scale_x_continuous(limits=c(1, 6), breaks=seq(1, 6)) +
    scale_y_continuous(limits=c(0, 155), breaks=seq(0, 150, 50), expand=c(0, 0)) +
    xlab("") +
    theme(axis.title.y=element_blank()),
  fig_cv[["mifs"]] +
    scale_y_continuous(limits=c(-1820, -1180), breaks=seq(-1800, -1200, 200)) +
    xlab(""),
  fig_cv[["dada"]] +
    scale_y_continuous(limits=c(-465, -395), breaks=seq(-460, -400, 20)) +
    theme(axis.title.y=element_blank()),
  fig_cv[["uno5"]] +
    scale_y_continuous(limits=c(-605, -475), breaks=seq(-600, -480, 40)) +
    xlab("") +
    theme(axis.title.y=element_blank()),
  fig_pdf2[["mifs"]] +
    scale_x_continuous(limits=c(1, 6), breaks=seq(1, 6)) +
    scale_y_continuous(limits=c(0, 510), breaks=seq(0, 500, 100), expand=c(0, 0)) +
    theme(legend.position=c(.96, .85)) +
    xlab(""),
  fig_pdf2[["dada"]] +
    scale_x_continuous(limits=c(1, 6), breaks=seq(1, 6)) +
    scale_y_continuous(limits=c(0, 155), breaks=seq(0, 150, 50), expand=c(0, 0)) +
    theme(axis.title.y=element_blank(),
          legend.position=c(.96, .85)),
  fig_pdf2[["uno5"]] + xlab("") +
    scale_x_continuous(limits=c(1, 6), breaks=seq(1, 6)) +
    scale_y_continuous(limits=c(0, 155), breaks=seq(0, 150, 50), expand=c(0, 0)) +
    theme(axis.title.y=element_blank(),
          legend.position=c(.96, .85)),
  align="hv", nrow=3,
  label_x=c(.51, .55, .5, .51, .55, .5, .51, .55, .5), label_y=.97,
  labels=c("(a) NON-DN", "(b) DADA2", "(c) UNOISE3",
           "(d) NON-DN", "(e) DADA2", "(f) UNOISE3",
           "(g) NON-DN", "(h) DADA2", "(i) UNOISE3")
))
save_plot(paste0(path_output, "/07-Fig_publ.svg"), fig_publ,
          base_asp=1, ncol=3, nrow=3)

# Save data
## Save R objects
path_saved_object <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_object, recursive=TRUE)
saveRDS(thresh_tab, paste0(path_saved_object, "/thresh_tab", ".obj"))
saveRDS(fig_histo, paste0(path_saved_object, "/fig_histo", ".obj"))
# Save workspace and session info
save.image(paste0(path_output, "/gmmdn.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/gmmdn.info"))
