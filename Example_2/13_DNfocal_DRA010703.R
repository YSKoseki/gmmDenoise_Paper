# 13_DNfocal_DRA010703.R
# Last updated on 2023.6.14 by YK
# An R script to generate graphical representations of the GMM-based denoising results in the focal species, Mugil cephalus
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
path_input1 <- "./12-DNsmry_DRA010703/01-Saved_object/phylo_pre.obj"
## GMM-inferred ASV read count cutoff threshold
path_input2 <- "./12-DNsmry_DRA010703/01-Saved_object/phylo_post.obj"
## Fasta of reference sequences
path_ref <- "./13-M_cephalus_ref.fa"
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
plot.tree(phylo_pre, "dada", root.node = 1484, xmax = 13.2, text.size = .3)

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
plot.tree(phylo_pre2, "dada", xmax = 13.2, text.size = .3)

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
plot.tree(phylo_pre3, "dada", xmax = 2.6, text.size = .3)

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
    refseq(x) %>% find.identical(ref["NC003182_NWP1"], .)
  }) # Returns "M. cephalus [0008]"
phylo_pre3 %>% lapply(
  function(x) {
    refseq(x) %>% find.identical(ref["KM368340_NWP2"], .)
  }) # Returns "M. cephalus [0062]"
## Check sequence lengths, with arranging ASVs in decreasing order
options("showHeadLines" = 60); options("showTailLines" = 60)
phylo_pre3 %>% lapply(
  function(x) {
    dna <- speedyseq::filter_tax_table(x, abbtag == "M. cephalus") %>% refseq()
    dna <- dna[order(Biostrings::width(dna), decreasing = TRUE), ]
    return(dna)
  })
## Replace the ID names of the ASVs identical to the reference sequences
phylo_pre4 <- phylo_pre3 %>% lapply(
  function(x) {
    taxtab <- tax_table(x)
    idseqnam1 <- refseq(x) %>%
      find.identical(ref["NC003182_NWP1"], .) %>%
      names()
    idseqnam2 <- refseq(x) %>%
      find.identical(ref["KM368340_NWP2"], .) %>%
      names()
    idseq1 <- taxtab[idseqnam1, ]
    idseq2 <- taxtab[idseqnam2, ]
    oldnam1 <- idseq1[, "glob"]
    oldnam2 <- idseq2[, "glob"]
    taxtab[idseqnam1, "glob"] <- paste0(oldnam1, ", NC003182")
    taxtab[idseqnam2, "glob"] <- paste0(oldnam2, ", KM368340")
    taxa_names(x) <- taxtab[, "glob"] %>% as.vector()
    return(x)
  }
)
## Test-plot the trees
plot.tree(phylo_pre4, "dada", subset.node = 1358, text.size = 3, xmax = .07)
plot.tree(phylo_pre4, "uno5", subset.node = 1562, text.size = 3, xmax = .04)
plot.tree(phylo_pre4, "uno2", subset.node = 784, text.size = 3, xmax = .06)
plot.tree(phylo_pre4, "mifs", subset.node = 7401, text.size = 1.5, xmax = .04)
## Arrange the trees for publication
focalnode <- list(dada = 1358, uno5 = 1562, uno2 = 784, mifs = 7401)
linesize <- list(dada = .5, uno5 = .5, uno2 = .5, mifs = .2)
labsize <- list(dada = 4, uno5 = 4, uno2 = 4, mifs = 1.5)
xmax2 <- list(dada = .083, uno5 = .048, uno2 = .069, mifs = .0394)
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
tree_focal[["uno2"]] <- tree_focal[["uno2"]] %>%
  flip(1, 8) %>% flip(2, 9) %>% flip(4, 7)
## Annotate the major clades
annode <- list(dada = c(30, 51), uno5 = c(26, 41), uno2 = c(6, 8), mifs = c(121, 230))
annof <- list(dada = .055, uno5 = .032, uno2 = .046, mifs = .01)
annosize <- list(dada = 4.5, uno5 = 4.5, uno2 = 4.5, mifs = 3.5)
annotof <- list(dada = -0.0006, uno5 = -0.0004, uno2 = -0.0006, mifs = -0.0004)
## Save the plots
tree_focal2 <- list(tree_focal, annode, annof, annosize, annotof, datnam) %>%
  pmap(function(p, n, of, s, tof, t) {
    p2 <- p +
      geom_cladelab(node = n[1], barcolor = "grey80", barsize = 9, offset = of,
                    label = "NWP1", fontsize = s, offset.text = tof,
                    hjust = "center", angle = 90, align = TRUE, textcolour = "white") +
      geom_cladelab(node = n[1], barcolor = alpha("grey80", 0), barsize = 9, offset = of,
                    label = "NWP1", fontsize = s, offset.text = tof,
                    hjust = "center", angle = 90, align = TRUE) +
      geom_cladelab(node = n[2], barcolor = "grey80", barsize = 9, offset = of,
                    label = "NWP2", fontsize = s, offset.text = tof,
                    hjust = "center", angle = 90, align = TRUE) +
      geom_cladelab(node = n[2], barcolor = alpha("grey80", 0), barsize = 9, offset = of,
                    label = "NWP2", fontsize = s, offset.text = tof,
                    hjust = "center", angle = 90, align = TRUE)
    save_plot(paste0(path_output, "/03-Focaltree_", t, ".svg"), p2,
              base_height = 7, base_asp = 1 / 1)
    return(p2)
  })


# Plot tree with heatmap
xmax3 <- list(dada = .01, uno5 = .01, uno2 = .01, mifs = .01)
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
