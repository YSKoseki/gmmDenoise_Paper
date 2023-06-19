# 13_DNfocal_DRA005106.R
# Last updated on 2023.6.14 by YK
# An R script to generate graphical representations of the GMM-based denoising results in the focal species, Liobagrus reinii
# R 4.1.2

# Packages required
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9018
library(ape); packageVersion("ape") # 5.6.2
library(treeio); packageVersion("treeio") # 1.18.1
library(ggtree); packageVersion("ggtree") # 3.2.1
library(cowplot); packageVersion("cowplot") # 1.1.1
library(svglite); packageVersion("svglite") # 2.1.0
library(aplot); packageVersion("aplot") # 0.1.6
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Paths
## List of phyloseq objects
path_input1 <- "./12-DNsmry_DRA005106/01-Saved_object/phylo_pre.obj"
## GMM-inferred ASV read count cutoff threshold
path_input2 <- "./12-DNsmry_DRA005106/01-Saved_object/phylo_post.obj"
## Fasta of reference sequences
path_ref <- "./13-L_reinii_ref.fa"
## Output directory
path_output <- "./13-DNfocal_DRA005106"
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
plot.tree(phylo_pre, "dada", root.node = 557, xmax = 1.8, text.size = .3)

# Replace trees with mammal-rooted trees
rootnode <- list(dada = 557, uno5 = 662, uno2 = 319, mifs = 2386)
phylo_pre2 <- phylo_pre
for (i in 1:length(phylo_pre2)) {
  rootedtree <- NULL
  rootedtree <- phy_tree(phylo_pre2[[i]]) %>% ape::root(node = rootnode[[i]])
  phy_tree(phylo_pre2[[i]]) <- rootedtree
}
rm(i)

# Check if phyloseq trees were replaced with the mammal-rooted ones
plot.tree(phylo_pre2, "dada", xmax = 1.8)

# Update phyloseq taxon tables for better annotation and omitting non-fish (mammal) clade
fsclass <- c("Actinopteri")
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
plot.tree(phylo_pre3, "dada", xmax = .8, text.size = .3)

# Plot global trees with denoised ASVs annotated by color
datnam <- phylo_pre %>% names() %>% as.list()
xmax <- list(dada = .8, uno5 = 1.1, uno2 = 1.1, mifs = 1.0)
tree_glob <- list(phylo_pre3, xmax, datnam) %>% 
  pmap(function(x, y, z) {
    p <- x %>%
      ggtree(size = .1) +
      geom_text(aes(label = node), size = .3, hjust = -.2) +
      geom_tiplab(aes(color = denoised), size = .3, hjust = -.1) +
      scale_color_manual(values = c("black", "firebrick")) +
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
    save_plot(paste0(path_output, "/02-Globaltree_", z, "_alp.svg"), p,
              base_height = 7, base_asp = 1 / 1)
    return(p)
  })


# Plot trees of the focal species 
## Read reference sequence fasta
ref <- Biostrings::readDNAStringSet(path_ref)
## Find ASVs identical to the reference sequences
find.identical<- function(DNAString, DNAStringset, type = "complete") {
  ref <- DNAString %>% as.character()
  seqs <- DNAStringset %>% as.character()
  if (type == "complete") {
    rownum <- match(ref, seqs)
  }
  if (type == "partial") {
    rownum <- stringr::str_detect(ref, seqs) %>% which()
  }
  z <- DNAStringset[rownum]
  return(z)
}
phylo_pre3 %>% lapply(
  function(x) {
    refseq(x) %>% find.identical(ref["LC468894_Clade1"], .)
  }) # Returns "L. reinii [0071]"
phylo_pre3 %>% lapply(
  function(x, y) {
    refseq(x) %>% find.identical(ref["LC146147_Clade2"], .)
  }) # Returns "L. reinii [0037]"
## Check sequence lengths, with arranging ASVs in decreasing order
options("showHeadLines" = 10); options("showTailLines" = 10)
phylo_pre3 %>% lapply(
  function(x) {
    dna <- speedyseq::filter_tax_table(x, abbtag == "L. reinii") %>% refseq()
    dna <- dna[order(Biostrings::width(dna), decreasing = TRUE), ]
    return(dna)
  })
## Replace the ID names of the ASVs identical to the reference sequences
phylo_pre4 <- phylo_pre3 %>% lapply(
  function(x) {
    taxtab <- tax_table(x)
    idseqnam1 <- refseq(x) %>%
      find.identical(ref["LC468894_Clade1"], .) %>%
      names()
    idseqnam2 <- refseq(x) %>%
      find.identical(ref["LC146147_Clade2"], .) %>%
      names()
    idseq1 <- taxtab[idseqnam1, ]
    idseq2 <- taxtab[idseqnam2, ]
    oldnam1 <- idseq1[, "glob"]
    oldnam2 <- idseq2[, "glob"]
    taxtab[idseqnam1, "glob"] <- paste0(oldnam1, ", LC468894")
    taxtab[idseqnam2, "glob"] <- paste0(oldnam2, ", LC146147")
    taxa_names(x) <- taxtab[, "glob"] %>% as.vector()
    return(x)
  }
)
## Test-plot the trees
plot.tree(phylo_pre4, "dada", subset.node = 432, text.size = 4, xmax = .1)
plot.tree(phylo_pre4, "uno5", subset.node = 526, text.size = 4, xmax = .1)
plot.tree(phylo_pre4, "uno2", subset.node = 256, text.size = 4, xmax = .1)
plot.tree(phylo_pre4, "mifs", subset.node = 1921, text.size = 4, xmax = .1)
## Arrange the trees for publication
focalnode <- list(dada = 432, uno5 = 526, uno2 = 256, mifs = 1921)
linesize <- list(dada = .7, uno5 = .7, uno2 = .7, mifs = .7)
labsize <- list(dada = 5.5, uno5 = 5.5, uno2 = 5.5, mifs = 5.5)
xmax2 <- list(dada = .084, uno5 = .113, uno2 = .094, mifs = .11)
tree_focal <- list(phylo_pre4, focalnode, linesize, labsize, xmax2) %>%
  pmap(function(x, y, z, s, t) {
    tree <- x %>% phy_tree() %>% treeio::tree_subset(node = y, levels_back = 0)
    taxdf <- x %>% tax_table() %>% as_tibble()
    p <- ggtree(tree, size = z) %<+% taxdf
    p2 <- p +
      geom_tiplab(aes(color = denoised), size = s, hjust = -.05, align = TRUE,
                  linesize = z / 2) +
      scale_color_manual(values = c("black", "firebrick")) +
      theme_tree2(legend.position = "none") +
      xlim_tree(t)
    return(p2)
  })
## Flip tree nodes for clade consistency with Nakagawa et al.'s (2016) tree
tree_focal[["dada"]] <- tree_focal[["dada"]] %>%
  flip(1, 3) %>% flip(2, 5) %>% flip(8, 10)
tree_focal[["uno5"]] <- tree_focal[["uno5"]] %>%
  flip(12, 15)
tree_focal[["uno2"]] <- tree_focal[["uno2"]] %>%
  flip(1, 3) %>% flip(2, 5)
tree_focal[["mifs"]] <- tree_focal[["mifs"]] %>%
  flip(2, 5) %>% flip(16, 18) %>% flip(1, 27)
## Annotate the major clades
annode <- list(dada = c(8, 10), uno5 = c(12, 15), uno2 = c(7, 9), mifs = c(16, 18))
annof <- list(dada = .085, uno5 = .11, uno2 = .095, mifs = .11)
annosize <- list(dada = 5.5, uno5 = 5.5, uno2 = 5.5, mifs = 5.5)
annotof <- list(dada = -0.0005, uno5 = -0.0007, uno2 = -0.0007, mifs = -0.0007)
## Save the plots
tree_focal2 <- list(tree_focal, annode, annof, annosize, annotof, datnam) %>%
  pmap(function(p, n, of, s, tof, t) {
    p2 <- p +
      geom_cladelab(node = n[1], barcolor = "grey80", barsize = 9, offset = of,
                    label = "Clade 1", fontsize = s, offset.text = tof,
                    hjust = "center", angle = 90, align = TRUE, textcolour = "white") +
      geom_cladelab(node = n[1], barcolor = alpha("grey80", 0), barsize = 9, offset = of,
                    label = "Clade 1", fontsize = s, offset.text = tof,
                    hjust = "center", angle = 90, align = TRUE) +
      geom_cladelab(node = n[2], barcolor = "grey80", barsize = 9, offset = of,
                    label = "Clade 2", fontsize = s, offset.text = tof,
                    hjust = "center", angle = 90, align = TRUE) +
      geom_cladelab(node = n[2], barcolor = alpha("grey80", 0), barsize = 9, offset = of,
                    label = "Clade 2", fontsize = s, offset.text = tof,
                    hjust = "center", angle = 90, align = TRUE)
    save_plot(paste0(path_output, "/03-Focaltree_", t, ".svg"), p2,
              base_height = 7, base_asp = 1 / 1)
    return(p2)
  })


# Plot tree with heatmap
xmax3 <- list(dada = .115, uno5 = .155, uno2 = .13, mifs = .15)
tree_focal_heat <- list(phylo_pre4, tree_focal2, datnam, xmax3) %>%
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
      geom_tile(aes(fill = logabund), color = "white", size = 1) +
      scale_fill_gradient(low = "#fae7e7", high = "black", name = "Abund") +
      scale_y_discrete(expand = c(0, 0)) +
      theme(axis.title.x = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1),
            axis.text.y = element_blank(),
            axis.line.x = element_line(size = .5),
            axis.line.y = element_line(size = .5),
            axis.ticks.x = element_line(size = .5),
            axis.ticks.y = element_line(size = .5),
            panel.background = element_rect(fill = "white"),
            legend.title = element_text(size = rel(1.3)),
            legend.text = element_text(size = rel(1.3))) +
      xlab("Sample") +
      ylab(NULL)
    comp.p <- p2 %>% aplot::insert_left(p, width = .8)
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
saveRDS(phylo_pre4, paste0(path_saved_object, "/phylo_pre4", ".obj"))
saveRDS(tree_glob, paste0(path_saved_object, "/tree_glob", ".obj"))
## Save workspace and session info
save.image(paste0(path_output, "/dnfocal.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/dnfocal.info"))
