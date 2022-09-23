# 11_GMMdn_DRA005106.R
# Last updated on 2022.7.3 by YK
# An R script to infer read count cut-off threshold for ASVs, based on Gaussian mixture modeling (GMM)
# R 4.1.2

# Packages required
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9018
library(mixtools); packageVersion("mixtools") # 1.2.0
library(cowplot); packageVersion("cowplot") # 1.1.1
library(wesanderson); packageVersion("wesanderson") # 0.3.6
library(tidyverse); packageVersion("tidyverse") # 1.3.1
library(gmmDenoise); packageVersion("gmmDenoise") # 0.1.0

# Paths
## Input: the list of phyloseq objects
path_input <- "./10-Phyloseq_DRA005106/01-Saved_object/phylo.obj"
## Output directory
path_output <- "./11-GMMdn_DRA005106"

# Load the phyloseq objects
phylo <- readRDS(path_input)
reads <- phylo %>%
  lapply(function(x) {x %>% otu_table() %>% colSums()})

# Step 0: Plot ASV read count histogram (in log scale)
fig_histo <- reads %>% 
  lapply(function(x) asvhist(x, scale="log", ylim=c(0, 200)))
theme_set(cowplot::theme_cowplot())
(fig_histo_grid <- fig_histo %>% 
  cowplot::plot_grid(plotlist=., 
                     labels=toupper(names(fig_histo)), 
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
(fig_cv_grid <- fig_cv %>% 
    cowplot::plot_grid(plotlist=., 
                       labels=toupper(names(fig_cv)), 
                       label_x=1, label_y=1, hjust=2, vjust=2))
save_plot(paste0(path_output, "/03-Fig_cv_grid.svg"), fig_cv_grid, ncol=2, nrow=2)

# Unused
# ## Approach 2: Sequential parametric bootstrap tests
# ### The bootstrapping
# set.seed(102)
# paraboot <- reads %>% 
#   lapply(function(x) {
#     gmmbs(log10(x), B=5000, epsilon=1e-2)
#   })
# ### Plot the bootstrap test results
# fig_paraboot <- paraboot %>% lapply(autoplot)
# ### Flatten the hierarchical structure of 'fig_paraboot'
# flatten.fig.paraboot <- function(fig.paraboot) {
#   max.cols <- fig.paraboot %>% sapply(length) %>% max()
#   flat.h <- map2(fig.paraboot, max.cols,
#                    function(x, y) {
#                      n.figs <- length(x)
#                      n.cols <- y
#                      n.diff <- n.cols-n.figs
#                      if (n.diff>0) 
#                        z <- rep(list(NULL), n.diff) %>% c(x, .)
#                      else 
#                        z <- x
#                    }) %>% 
#     unlist(recursive=FALSE, use.names=FALSE)
#   # Give name labels of 'fig.paraboot' to the element of left-most plotted histograms
#   elms <- seq(1, length(flat.h), by=max.cols)
#   nams <- names(fig.paraboot)
#   names(flat.h)[elms] <- nams
#   return(flat.h)
# }
# fig_paraboot2 <- flatten.fig.paraboot(fig_paraboot)
# ### Display the plots in grid
# theme_set(cowplot::theme_cowplot())
# (fig_paraboot_grid <- fig_paraboot2 %>% 
#     cowplot::plot_grid(plotlist=.,
#                        ncol=length(fig_paraboot2)/length(fig_paraboot),
#                        nrow=length(fig_paraboot),
#                        labels=toupper(names(fig_paraboot2)), label_size=16,
#                        label_x=0, label_y=1, hjust=-.3, vjust=1.6))
# save_plot(paste0(path_output, "/04-Fig_paraboot_grid.svg"), 
#           fig_paraboot_grid, 
#           ncol=length(fig_paraboot2)/length(fig_paraboot),
#           nrow=length(fig_paraboot))
# ### Summary statistics of the bootstrap tests
# (summary_parboot <- paraboot %>% lapply(summary))

# Step 2: Model fitting with the selected number of components
## The selected number of components
k <- list(data=2, uno5=3, uno2=1, mifs=3)
## Fit k-component GMMs
set.seed(103)
gmm_fit <- map2(reads, k, 
            function(x, k) gmmem(log10(x), k, maxit=5000, epsilon=1e-2))

# Step 3: Plot GMM-estimated count density functions (CDF) on histogram
## Draw the plots
fig_pdf <- gmm_fit %>%
  lapply(function(x) {
       autoplot(x, xlim=c(1, 6), ylim=c(0, 200))
     })
theme_set(cowplot::theme_cowplot())
fig_pdf_grid <- fig_pdf %>% 
  cowplot::plot_grid(plotlist=., 
                     labels=toupper(names(fig_pdf)), 
                     label_x=.75, label_y=1, hjust=0, vjust=2)
save_plot(paste0(path_output, "/05-Fig_pdf_grid.svg"), fig_pdf_grid, ncol=2, nrow=2)

# Step 4: Infer read count threshold for ASV cut-off
## Calculate 95-percentiles of individual mixture components of fitted GMMs
lower95_log <- gmm_fit %>% 
  lapply(function(x) {
    ncomp <- x$mu %>% length()
    quantile(x, comp=1:ncomp)
  })
(thresh_tab <- list(dada=c(1, NA), uno5=c(NA, 1, NA), uno2=NA, mifs=c(NA, 1, NA)) %>% 
    map2(lower95_log, function(x, y) x*y) %>% 
    lapply(function(x) {
      if (all(is.na(x))) NA
      else na.omit(x)
    }) %>% 
    unlist() %>% 
    data.frame(log=.) %>% 
    rownames_to_column(var="data") %>% 
    mutate(norm=ceiling(10^log)))
write.csv(thresh_tab, paste0(path_output, "/06-thresh_tab.csv"))
## Draw vertical lines indicating 95-percentiles to the plots
fig_pdf2 <- list(dada=c(1, NA), uno5=c(NA, 1, NA), uno2=NA, mifs=c(NA, 1, NA)) %>% 
  map2(lower95_log, function(x, y) x*y) %>% 
  map2(gmm_fit, .,
       function(x, y) {
         autoplot(x, vline=y, xlim=c(1, 6), ylim=c(0, 200))
       })
theme_set(cowplot::theme_cowplot())
fig_pdf_grid2 <- fig_pdf2 %>% 
    cowplot::plot_grid(plotlist=., 
                       labels=toupper(names(fig_pdf)), 
                       label_x=.75, label_y=1, hjust=0, vjust=2)
save_plot(paste0(path_output, "/07-Fig_pdf_grid2.svg"), fig_pdf_grid2, ncol=2, nrow=2)

# Save data
## Save R objects
path_saved_object <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_object, recursive=TRUE)
saveRDS(thresh_tab, paste0(path_saved_object, "/thresh_tab", ".obj"))
# Save workspace and session info
save.image(paste0(path_output, "/gmmdn.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/gmmdn.info"))