# 09_GMMdn_DRA006638.R
# Last updated on 2025.6.7 by YK
# An R script to infer true ASVs by the denoising method based on Gaussian mixture modeling (GMM)
# R 4.1.2

# Packages required
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9018
library(cowplot); packageVersion("cowplot") # 1.1.1
library(wesanderson); packageVersion("wesanderson") # 0.3.6
library(tidyverse); packageVersion("tidyverse") # 1.3.1
library(gmmDenoise); packageVersion("gmmDenoise") # 0.3.1

# Paths
## Input: the list of phyloseq objects
path_input <- "./08-Phyloseq_DRA006638/01-Saved_object/phylo.obj"
## Output directory
path_output <- "./09-GMMdn_DRA006638"
dir.create(path_output, recursive=TRUE)

# Load the phyloseq objects and generate read data
phylo <- readRDS(path_input)
reads <- phylo %>%
  lapply(function(x) {x %>% otu_table() %>% colSums()})
nrepl <- phylo %>%
  lapply(function(x) {
    y <- x %>% otu_table()
    y[y>0] <- 1
    z <- y %>% colSums()
    return(z)})
truehap <- phylo %>%
  lapply(function(x) {
    y <- x %>% tax_table() %>% as.data.frame()
    y$istruehap <- factor(y$istruehap, levels=c("True", "False"))
    z <- y  %>% rownames_to_column(var="asv") %>%
      select(asv, truehap, istruehap) %>% 
      as_tibble()
    return(z)})
reads_tib <- list(reads, nrepl, truehap) %>%
  pmap(function(x, y, z) {
    nam <- names(x)
    v <- tibble(asv=nam, reads=x, nrepl=y)
    w <- v %>% left_join(z, by="asv")
    return(w)})
(truefalse_pre <- reads_tib %>%
  lapply(function(x) {
    x %>% group_by(istruehap) %>% summarize(n=n())
  }))
(detect_in_all <- reads_tib %>%
  lapply(function(x) {
    x %>% group_by(istruehap) %>%
      filter(nrepl==15) %>%
      summarize(detect_in_all=n_distinct(asv))
  }))
write.csv(c(truefalse_pre[c(3, 1, 2)], detect_in_all[c(3, 1, 2)]),
          paste0(path_output, "/02-truefalse_pre.csv"))

# Plot of reads vs. detection rates (as Figure 3 in Tsuji et al. 2020)
fig_scatter <- reads_tib %>%
  lapply(function(x) {
    y <- x %>%
      ggplot(aes(x=nrepl, y=log10(reads), color=istruehap)) +
      geom_point(shape=21, size=3) +
      scale_color_manual(values=c("black", "firebrick")) +
      scale_x_continuous(limits=c(1, 15), breaks=seq(1, 15, 2)) +
      scale_y_continuous(limits=c(0, 6), breaks=seq(0, 6, 1)) +
      xlab("Detection rate in 15 PCR replicates") +
      ylab("log10(Total reads)") +
      guides(color="none")
  })
theme_set(cowplot::theme_cowplot())
(fig_scatter_grid <- cowplot::plot_grid(
  fig_scatter[["nodn"]] + xlab(""),
  fig_scatter[["dada"]] + theme(axis.title.y=element_blank()),
  fig_scatter[["uno2"]] + xlab("") + theme(axis.title.y=element_blank()),
  nrow=1, label_x=c(.13, .07, .06), label_y=.95, hjust=-.2,
  labels=c("(a) NON-DN", "(b) DADA2", "(c) UNOISE3")
))
save_plot(paste0(path_output, "/03-Fig_scatter_grid.svg"), fig_scatter_grid,
          base_asp=0.9, ncol=3, nrow=1)

# Histogram plots
fig_histo <- reads %>% 
  lapply(function(x) asvhist(x, scale="log", type="freq", nbins="Sturges",
                             xlim = c(0, 6)))
theme_set(cowplot::theme_cowplot())
(fig_histo_grid <- cowplot::plot_grid(
  fig_histo[["nodn"]] +
    scale_y_continuous(limits=c(0, 800), expand=c(0, 0)) +
    xlab(""),
  fig_histo[["dada"]] +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    theme(axis.title.y=element_blank()),
  fig_histo[["uno2"]] +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    xlab("") +
    theme(axis.title.y=element_blank()),
  nrow=1, label_x=c(.54, .55, .52), label_y=.97, hjust=-.35,
  #nrow=1, label_x=c(.52, .59, .48), label_y=.98, hjust=-.80,
  labels=c("NON-DN", "DADA2", "UNOISE3")
))
save_plot(paste0(path_output, "/04-Fig_histo_grid.svg"), fig_histo_grid,
          base_asp=0.9, ncol=3, nrow=1)

# Cross-validation analysis for selecting the number of components
set.seed(006638)
crossval <- reads %>% 
  lapply(function(x) gmmcv(log10(x), maxk=5, maxit=5000, maxrestarts=1000))
fig_cv <- crossval %>% lapply(autoplot, type="density")
theme_set(cowplot::theme_cowplot())
(fig_cv_grid <- cowplot::plot_grid(
  fig_cv[["nodn"]] + xlab(""),
  fig_cv[["dada"]] + theme(axis.title.y=element_blank()),
  fig_cv[["uno2"]] + xlab("") + theme(axis.title.y=element_blank()),
  nrow=1, label_x=c(.08, .59, .48), label_y=.98, hjust=-.80,
  labels=c("NON-DN", "DADA2", "UNOISE3")
))
save_plot(paste0(path_output, "/05-Fig_cv_grid.svg"), fig_cv_grid,
          base_asp=0.9, ncol=3, nrow=1)

# Fit k-component GMMs
k <- list(data=2, uno2=2, nodn=4)
gmm_fit <- map2(reads, k, 
            function(x, y) gmmem(log10(x), k=y, maxit=5000, maxrestarts=1000))

# Plot GMM-estimated count density functions (CDF) on histogram
fig_pdf <- gmm_fit %>%
  lapply(function(x) {
    autoplot(x, type="freq", nbins="Sturges", xlim=c(0, 6))
  })
theme_set(cowplot::theme_cowplot())
(fig_pdf_grid <- cowplot::plot_grid(
  fig_pdf[["nodn"]] + xlab("") +
    scale_y_continuous(limits=c(0, 800), expand=c(0, 0)) +
    theme(legend.position=c(.9, .9)) +
    xlab(""),
  fig_pdf[["dada"]] +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    theme(axis.title.y=element_blank(),
          legend.position=c(.9, .9)),
  fig_pdf[["uno2"]] +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    xlab("") +
    theme(axis.title.y=element_blank(),
          legend.position=c(.9, .9)),
  nrow=1, label_x=c(.54, .55, .52), label_y=.97, hjust=-.35,
  labels=c("NON-DN", "DADA2", "UNOISE3")
))
save_plot(paste0(path_output, "/06-Fig_pdf_grid.svg"), fig_pdf_grid,
          base_asp=0.9, ncol=3, nrow=1)

# Infer read count threshold for ASV cut-off
lower95_log <- gmm_fit %>% 
  lapply(function(x) {
    ncomp <- x$mu %>% length()
    quantile(x)
  })
thresh_tab <- lower95_log %>% 
    unlist() %>% 
    data.frame(log=.) %>% 
    rownames_to_column(var="data") %>% 
    mutate(norm=ceiling(10^log))
write.csv(thresh_tab[c(3, 1, 2), ], paste0(path_output, "/07-thresh_tab.csv"))

## Draw vertical lines indicating 95-percentiles to the pdf plots
fig_pdf2 <- lower95_log %>% 
  map2(gmm_fit, .,
       function(x, y) {
         autoplot(x, type="freq", nbins="Sturges", vline=y, xlim=c(0, 6))
       })
theme_set(cowplot::theme_cowplot())
fig_pdf_grid2 <- cowplot::plot_grid(
  fig_pdf2[["nodn"]] +
    scale_y_continuous(limits=c(0, 800), expand=c(0, 0)) +
    theme(legend.position=c(.94, .9)) +
    xlab(""),
  fig_pdf2[["dada"]] +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    theme(axis.title.y=element_blank(),
          legend.position=c(.43, .9)),
  fig_pdf2[["uno2"]] + xlab("") +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    theme(axis.title.y=element_blank(),
          legend.position=c(.94, .9)),
  nrow=1, label_x=c(.57, .14, .56), label_y=.97, hjust=-.35,
  labels=c("NON-DN", "DADA2", "UNOISE3")
)
save_plot(paste0(path_output, "/08-Fig_pdf_grid2.svg"), fig_pdf_grid2,
          base_asp=0.9, ncol=3, nrow=1)

# A summary table of the denoising effect
hapgroup_lev <- c("True pos", "True neg",
                  "False pos", "False neg")
reads_tib <- thresh_tab[, "norm"] %>%
  as.list() %>%
  map2(reads_tib, ., function(x, y) {
    z <- x %>%
      mutate(hapgroup=case_when(reads>y & istruehap=="True" ~
                                    factor("True pos", levels=hapgroup_lev),
                                  reads<=y & istruehap=="True" ~
                                    factor("False neg", levels=hapgroup_lev),
                                  reads>y & istruehap=="False" ~
                                    factor("False pos", levels=hapgroup_lev),
                                  reads<=y & istruehap=="False" ~
                                    factor("True neg", levels=hapgroup_lev)))
    return(z)
  })
(truefalse_post <- reads_tib %>%
  lapply(function(x) {
    x %>%
      group_by(hapgroup, .drop=FALSE) %>%
      summarize(n=n())
  }))
truefalse_post[c(3, 1, 2)] %>% as_tibble()
write.csv(truefalse_post[c(3, 1, 2)], paste0(path_output, "/09-truefalse_post.csv"))

# Visual representation of the denoising effect
fig_scatter2 <- thresh_tab[, "log"] %>%
  as.list() %>%
  map2(reads_tib, ., function(x, y) {
    z <- x %>%
      ggplot(aes(x=nrepl, y=log10(reads), color=hapgroup)) +
      annotate("rect", xmin=-Inf, xmax=14.5, ymin=-Inf, ymax=Inf,
               fill="black", alpha=.1) +
      geom_vline(xintercept=14.5, size=.05, color="grey80") +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=y,
               fill="firebrick", alpha=.1) +
      geom_hline(yintercept=y, size=.05, color="firebrick") +
      geom_point(shape=21, size=3) +
      scale_color_manual(values=c("black", "firebrick", "firebrick")) +
      scale_x_continuous(limits=c(1, 15), breaks=seq(1, 15, 2)) +
      scale_y_continuous(limits=c(0, 6), breaks=seq(0, 6, 1)) +
      xlab("Detection rate in 15 PCR replicates") +
      ylab("log10(Total reads)") +
      guides(color="none")
  })
theme_set(cowplot::theme_cowplot())
(fig_scatter_grid2 <- cowplot::plot_grid(
  fig_scatter2[["nodn"]] + xlab(""),
  fig_scatter2[["dada"]] + theme(axis.title.y=element_blank()),
  fig_scatter2[["uno2"]] + xlab("") + theme(axis.title.y=element_blank()),
  nrow=1, label_x=c(.13, .07, .06), label_y=.95, hjust=-.2,
  labels=c("(a) NON-DN", "(b) DADA2", "(c) UNOISE3")
))
save_plot(paste0(path_output, "/10-Fig_scatter_grid2.svg"), fig_scatter_grid2,
          base_asp=0.9, ncol=3, nrow=1)

# Compound plot for publication
theme_set(cowplot::theme_cowplot())
(fig_publ <- cowplot::plot_grid(
  fig_histo[["nodn"]] +
    scale_y_continuous(limits=c(0, 800), expand=c(0, 0)) +
    xlab(""),
  fig_histo[["dada"]] +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    theme(axis.title.y=element_blank()),
  fig_histo[["uno2"]] +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    xlab("") +
    theme(axis.title.y=element_blank()),
  fig_cv[["nodn"]] +
    scale_y_continuous(limits=c(-600, -200)) +
    xlab("") +
    annotate("text", x=4, y=-238, label="*", size=10),
  fig_cv[["dada"]] +
    scale_y_continuous(limits=c(-44, -37)) +
    theme(axis.title.y=element_blank()) +
    annotate("text", x=2, y=-38, label="*", size=10),
  fig_cv[["uno2"]] +
    scale_y_continuous(limits=c(-150, 0)) +
    xlab("") +
    theme(axis.title.y=element_blank()) +
    annotate("text", x=2, y=-17, label="*", size=10),
  fig_pdf2[["nodn"]] +
    scale_y_continuous(limits=c(0, 800), expand=c(0, 0)) +
    theme(legend.position=c(.94, .9)) +
    xlab(""),
  fig_pdf2[["dada"]] +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    theme(axis.title.y=element_blank(),
          legend.position=c(.37, .9)),
  fig_pdf2[["uno2"]] + xlab("") +
    scale_y_continuous(limits=c(0, 20), expand=c(0, 0)) +
    theme(axis.title.y=element_blank(),
          legend.position=c(.41, .9)),
  align="hv", nrow=3,
  label_x=c(.089, .098, .083, .089, .098, .083, .5, .098, .083), label_y=.97,
  labels=c("(a) NON-DN", "(b) DADA2", "(c) UNOISE3",
           "(d) NON-DN", "(e) DADA2", "(f) UNOISE3",
           "(g) NON-DN", "(h) DADA2", "(i) UNOISE3")
))
save_plot(paste0(path_output, "/11-Fig_publ.svg"), fig_publ,
          base_asp=1, ncol=3, nrow=3)

# Save data
## Save R objects
path_saved_object <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_object, recursive=TRUE)
saveRDS(reads_tib, paste0(path_saved_object, "/reads_tib", ".obj"))
saveRDS(thresh_tab, paste0(path_saved_object, "/thresh_tab", ".obj"))
# Save workspace and session info
save.image(paste0(path_output, "/gmmdn.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/gmmdn.info"))
