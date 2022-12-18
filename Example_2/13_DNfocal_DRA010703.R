# 13_DNfocal_DRA005106.R
# Last updated on 2022.12.18 by YK
# An R script to generate graphical representations of the GMM-based denoising results in a focal species
# R 4.1.2

# Packages required
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9018
library(cowplot); packageVersion("cowplot") # 1.1.1
library(svglite); packageVersion("svglite") # 2.1.0
library(ape); packageVersion("ape") # 5.6.2
library(treeio); packageVersion("treeio") # 1.18.1
library(ggtree); packageVersion("ggtree") # 3.2.1
library(ggmsa); packageVersion("ggmsa") # 1.3.3
library(aplot); packageVersion("aplot") # 0.1.6
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Paths
## List of phyloseq objects
path_input1 <- "./12-DNsmry_DRA010703/01-Saved_object/phylo_pre.obj"
## GMM-inferred ASV read count cutoff threshold
path_input2 <- "./12-DNsmry_DRA010703/01-Saved_object/phylo_post.obj"
## Output directory
path_output <- "./13-DNfocal_DRA010703"
dir.create(path_output, recursive = TRUE)

# Load the phyloseq objects
phylo_pre <- readRDS(path_input1)
phylo_post <- readRDS(path_input2)

# Test plotting of phyloseq trees
plot.tree <- function(phyloseq.list, list.id,
                      root.node = NULL, subset.node = NULL,
                      text.size = .3, xmax = 1.0) { 
  tree <- phyloseq.list[[list.id]] %>% phy_tree()
  if(!is.null(root.node)) {
    tree <- tree %>% ape::root(node = root.node)
  }
  if(!is.null(subset.node)) {
    tree <- tree %>% treeio::tree_subset(node = subset.node, levels_back = 0)
  }
  ggtree(tree, size = .1) +
    geom_text(aes(label = node), size = text.size, hjust = -.2) +
    geom_tiplab(size = text.size, hjust = -.2) +
    theme_tree2() +
    xlim(0, xmax)
}
plot.tree(phylo_pre, "mifs", root.node = NULL, subset.node = NULL,
          xmax = 2, text.size = .2)

# Replace trees with mammal-rooted trees
rootnode <- list(dada = 1484, uno5 = 1709, uno2 = 506, mifs = 7689)
phylo_pre2 <- phylo_pre
for (i in 1:length(phylo_pre2)) {
  rootedtree <- NULL
  rootedtree <- phy_tree(phylo_pre2[[i]]) %>% ape::root(node = rootnode[[i]])
  phy_tree(phylo_pre2[[i]]) <- rootedtree
}
rm(i)

# Check if phyloseq trees were replaced with the mammal-rooted ones
plot.tree(phylo_pre2, "mifs", xmax = 2, text.size = .2)

# Update phyloseq taxon tables for better annotation and omitting non-fish (mammal) clade
fsclass <- c("Actinopteri", "Chondrichthyes")
phylo_pre3 <- phylo_pre2 %>%
  map2(phylo_post,
       function(x, y) {
         # Identify denoised ASVs
         asv.x <- x %>% taxa_names()
         asv.y <- y %>% taxa_names()
         asv.deno <- setdiff(asv.x, asv.y)
         # Add a column for identification of denoised ASVs ("yes" or "no")
         new.x <- x %>%
           speedyseq::mutate_tax_table(
             denoised = if_else(.otu %in% asv.deno,
                                factor("yes", levels = c("yes", "no")),
                                factor("no", levels = c("yes", "no")))) %>%
           speedyseq::select_tax_table(.otu, denoised, everything()) %>%
           speedyseq::filter_tax_table(class %in% fsclass)
         return(new.x)
       })

# Check if non-fish (mammal) clades were omitted from phyloseq trees
plot.tree(phylo_pre3, "mifs", text.size = .2, xmax = 1.1)

# Plot global trees with denoised ASVs annotated by color
datnam <- phylo_pre %>% names() %>% as.list()
xmax <- list(dada = 2.6, uno5 = 1.7, uno2 = 2.5, mifs = 1.1)
tree_glob <- list(phylo_pre3, xmax, datnam) %>% 
  pmap(function(x, y, z) {
    p <- x %>%
      ggtree(size = .1) +
      geom_text(aes(label = node), size = .3, hjust = -.2) +
      geom_tiplab(aes(color = denoised), size = .3, hjust = -.1) +
      scale_color_manual(values = c("black", "firebrick")) +
      scale_alpha_manual(values = c(1, 0)) +
      theme_tree(legend.position = "none") +
      xlim_tree(y)
    save_plot(paste0(path_output, "/02-Globaltree_", z, ".svg"), p,
              base_height = 7, base_asp = 1 / 1)
    return(p)
  })

# Another global tree plots, where denoised ASVs are not annotated for readability 
tree_glob_alp <- list(phylo_pre3, xmax, datnam) %>% 
  pmap(function(x, y, z) {
    p <- x %>%
      ggtree(size = .1) +
      geom_text(aes(label = node), size = .3, hjust = -.2) +
      geom_tiplab(aes(alpha = denoised), size = .3, hjust = -.1) +
      scale_alpha_manual(values = c(1, 0)) +
      theme_tree(legend.position = "none") +
      xlim_tree(y)
    #save_plot(paste0(path_output, "/02-Globaltree_", z, ".svg"), p,
    save_plot(paste0(path_output, "/02-Globaltree_", z, "_alp.svg"), p,
              base_height = 7, base_asp = 1 / 1)
    return(p)
  })

# Plot trees of a focal species 
## Test plot for finding the node of focal clade
phylo_pre3 <- phylo_pre3 %>% lapply(
  function(x) {
    taxa_names(x) <- tax_table(x) %>% `[`(, "glob") %>% as.vector()
    return(x)
  }
)
plot.tree(phylo_pre3, "mifs", subset.node = 7401, text.size = 2, xmax = .05)

## Publication-quality tree plot 
focalnode <- list(dada = 1358, uno5 = 1562, uno2 = 784, mifs = 7401)

linesize <- list(dada = .5, uno5 = .5, uno2 = .5, mifs = .2)
labsize <- list(dada = 4, uno5 = 4, uno2 = 5, mifs = 1.5)
xmax2 <- list(dada = .07, uno5 = .0407, uno2 = .0608, mifs = .037)
tree_focal <- list(phylo_pre3, focalnode, linesize, labsize, xmax2, datnam) %>%
  pmap(function(x, y, z, s, t, u) {
    tree <- x %>% phy_tree() %>% treeio::tree_subset(node = y, levels_back = 0)
    taxdf <- x %>% tax_table() %>% as_tibble()
    p <- ggtree(tree, size = z) %<+% taxdf
    p2 <- p +
      geom_tiplab(aes(color = denoised), size = s, hjust = -.05, align = TRUE,
                  linesize = z / 2) +
      scale_color_manual(values = c("black", "firebrick")) +
      theme_tree2(legend.position = "none") +
      xlim_tree(t)
    save_plot(paste0(path_output, "/03-Focaltree_", u, ".svg"), p2,
              base_height = 7, base_asp = 1 / 1)
    return(p2)
  })

## Flip the uno2 tree nodes for consistency
tree_focal[["uno2"]] + geom_text(aes(label = node))
tree_focal[["uno2"]] <- tree_focal[["uno2"]] %>% flip(2, 9)
tree_focal[["uno2"]] <- tree_focal[["uno2"]] %>% flip(1, 8)
tree_focal[["uno2"]] <- tree_focal[["uno2"]] %>% flip(4, 7)
save_plot(paste0(path_output, "/03-Focaltree_uno2.svg"),
          tree_focal[["uno2"]],
          base_height = 7, base_asp = 1 / 1)

# Plot tree with heatmap
xmax3 <- list(dada = .078, uno5 = .045, uno2 = .072, mifs = .0385)
tree_focal_heat <- list(phylo_pre3, tree_focal, datnam, xmax3) %>%
  pmap(function(x, y, z, s) {
    focal <- y %>% `[[`("data") %>% filter(isTip == TRUE) %>% select(label) %>%
      pull()
    df <- x %>% otu_table() %>% as_tibble() %>%
      mutate(label = .otu, sample = .sample, abund = .abundance,
             logabund = log10(.abundance + 1)) %>%
      select(label, sample, abund, logabund) %>% filter(label %in% focal)
    goodsample <- df %>% group_by(sample) %>% summarize(totab = sum(abund)) %>%
      filter(totab > 0) %>% select(sample) %>% pull()
    df <- df %>% filter(sample %in% goodsample)
    p <- y + xlim_tree(s) + theme_tree(legend.position = "none")
    p2 <- ggplot(df, aes(x = sample, y = label)) +
      geom_tile(aes(fill = logabund), color = "white", size = .5) +
      scale_fill_gradient(low = "#fae7e7", high = "black", name = "Abund") +
      scale_y_discrete(expand = c(0, 0)) +
      theme(axis.title.x = element_text(size = rel(1.2)),
            axis.text.x = element_text(size = rel(1.2), angle = 45, hjust = 1),
            axis.text.y = element_blank(),
            axis.line.x = element_line(size = .5),
            axis.line.y = element_line(size = .5),
            axis.ticks.x = element_line(size = .5),
            axis.ticks.y = element_line(size = .5),
            panel.background = element_rect(fill = "white"),
            legend.title = element_text(size = rel(1.1)),
            legend.text = element_text(size = rel(1.1))) +
      xlab("Sample") +
      ylab(NULL)
    comp.p <- p2 %>% aplot::insert_left(p, width = .7)
    save_plot(paste0(path_output, "/04-FocaltreeHeat_", z, ".svg"), comp.p,
              base_height = 7, base_asp = 1.618 / 1)
    return(comp.p)
  })

# Summary stats of haplotypes
nhap <- map2(tree_focal, datnam,
       function(x, y) {
         tib <- x %>% `[[`("data") %>%
           filter(isTip) %>% filter(!is.na((!!y))) %>%
           select(denoised, (!!y))
         tot <- tib %>% nrow()
         ret <- tib %>% filter(denoised == "no") %>% nrow()
         df <- data.frame(data = y, 
                          pre = tot, post = ret, retain = ret / tot)
       })
nhap %>% bind_rows() %>%
  write.csv(paste0(path_output, "/05-dn_stat.csv"))

# Save data
## Save R objects
path_saved_object <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_object, recursive = TRUE)
saveRDS(focalnode, paste0(path_saved_object, "/focalnode", ".obj"))
saveRDS(phylo_pre3, paste0(path_saved_object, "/phylo_pre3", ".obj"))
saveRDS(tree_glob, paste0(path_saved_object, "/tree_glob", ".obj"))
# Save workspace and session info
save.image(paste0(path_output, "/dnfocal.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/dnfocal.info"))
