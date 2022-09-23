# 13_DNfocal_DRA010703.R
# Last updated on 2022.7.5 by YK
# An R script to generate graphical representations of the GMM-based denoising results in a focal species
# R 4.1.2

# Packages required
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9018
library(cowplot); packageVersion("cowplot") # 1.1.1
library(ape); packageVersion("ape") # 5.6.2
library(treeio); packageVersion("treeio") # 1.18.1
library(ggtree); packageVersion("ggtree") # 3.2.1
library(ggmsa); packageVersion("ggmsa") # 1.3.3
library(aplot); packageVersion("aplot") # 0.1.6
library(Biostrings); packageVersion("Biostrings") # 2.62.0
library(pegas); packageVersion("pegas") # 1.1
library(svglite); packageVersion("svglite") # 2.1.0
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

# Update phyloseq taxon tables for better annotation and omitting mammal clade
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

# Check if non-fish clades were omitted from phyloseq trees
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
plot.tree(phylo_pre3, "mifs", subset.node = 7401, text.size = 2, xmax = .05)

## Publication-quality tree plot 
focalnode <- list(dada = 1358, uno5 = 1562, uno2 = 784, mifs = 7401)
linesize <- list(dada = .5, uno5 = .5, uno2 = .5, mifs = .2)
labsize <- list(dada = 4, uno5 = 4, uno2 = 5, mifs = 1.5)
xmax2 <- list(dada = .082, uno5 = .047, uno2 = .073, mifs = .04)
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

# Plot tree with multiple sequence alignment plot
xmax3 <- list(dada = .105, uno5 = .061, uno2 = .11, mifs = .042)
tree_focal_msa <- list(phylo_pre3, tree_focal, datnam, xmax3) %>%
  pmap(function(x, y, z, s) {
    focal <- y %>% `[[`("data") %>% filter(isTip == TRUE) %>% select(label) %>%
      pull()
    dnaseq <- x %>% refseq() %>% .[names(.) %in% focal] %>% 
      DECIPHER::AlignSeqs(processors = NULL) %>%
      DECIPHER::StaggerAlignment(processors = NULL) %>%
      ggmsa::tidy_msa()
    p <- y +
      xlim_tree(s) +
      geom_facet(panel = "MSA", data = dnaseq, geom = geom_msa,
                 font = NULL, color = "Chemistry_NT", border = "white") +
      scale_x_continuous(expand = c(.003, .003)) +
      theme_tree2(legend.position = "none",
                  strip.text = element_blank(),
                  panel.spacing = unit(0.2, 'cm'))
    p2 <- facet_widths(p, widths = c(.5, 1))
    save_plot(paste0(path_output, "/04-FocaltreeMSA_", z, ".svg"), p2,
              base_height = 7, base_asp = 1.618 / 1)
    return(p)
  })

# Plot tree with heatmap
xmax4 <- list(dada = .115, uno5 = .066, uno2 = .125, mifs = .043)
tree_focal_heat <- list(phylo_pre3, tree_focal, datnam, xmax4) %>%
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
      geom_tile(aes(fill = logabund), color = "white", size = .2) +
      scale_fill_gradient(low = "#fae7e7", high = "black", name = "Abund") +
      theme(axis.text.x = element_text(size = rel(.5), angle = 90, hjust = 1),
            axis.text.y = element_blank(),
            axis.line.x = element_line(size = .2),
            axis.line.y = element_line(size = .2),
            axis.ticks.x = element_line(size = .2),
            axis.ticks.y = element_line(size = .2)) +
      xlab("Sample (site)") +
      ylab(NULL)
    comp.p <- p2 %>% aplot::insert_left(p, width = .5)
    save_plot(paste0(path_output, "/05-FocaltreeHeat_", z, ".svg"), comp.p,
              base_height = 7, base_asp = 1.618 / 1)
    return(comp.p)
  })

# Haplotype network plots
## Create objects for plotting
hapnet <- phylo_pre3 %>%
  map2(focalnode, 
       function(x, y) {
         label <- x %>% phy_tree() %>%
           treeio::tree_subset(node = y, levels_back = 0) %>%
           `[[`("tip.label")
         phyfocal <- x %>% prune_taxa(label, .)
         taxa_names(phyfocal) <- phyfocal %>% tax_table() %>%
           `[`(, "glob") %>% as.vector()
         nsamp <- phyfocal %>% otu_table() %>%
           apply(2, function(x) {which(x > 0) %>% length()}) 
         freq <- phyfocal %>% otu_table() %>% colSums()
         dnabin <- phyfocal %>% refseq() %>%
           DECIPHER::AlignSeqs(processors = NULL, verbose = FALSE) %>%
           DECIPHER::StaggerAlignment(processors = NULL, verbose = FALSE) %>%
           rep(times = freq) %>% as.DNAbin()
         hap <- pegas::haplotype(dnabin, strict = TRUE,
                                 labels = names(freq))
         net <- pegas::haploNet(hap)
         dn_tib <- phyfocal %>% tax_table() %>%
           as_tibble() %>% select(.otu, denoised) %>%
           mutate(color = if_else(denoised == "no", 1, 2)) %>%
           select(.otu, color)
         dn <- dn_tib %>% select(color) %>% pull()
         names(dn) <- dn_tib %>% select(.otu) %>% pull()
         # See if I haven't made any mistakes
         freq_net <- attr(net, "freq") %>% as.numeric()
         names(freq_net) <- attr(net, "labels")
         if (identical(freq, freq_net)) {
           cat("The ASV frequency data are arranged consistent with the ASV names.\n")
         } else {
           stop("The ASV frequency data are NOT arranged consistent with the ASV names!\n")
         }
         if (identical(names(nsamp), names(freq_net))) {
           cat("The ASV-positive sample size data are arranged consistent with the ASV names.\n")
         } else {
           stop("The ASV-positive sample size data are NOT arranged consistent with the ASV names!\n")
         }
         if (identical(names(dn), names(freq_net))) {
           cat("The denoising identities of ASVs are arranged consistent with the ASV names.\n")
         } else {
           stop("The denoising identities of ASVs are NOT arranged consistent with the ASV names!\n")
         }
         z <- list(net = net, nsamp = nsamp, dn = dn)
         return(z)
       })
## Set XY coordinates of haplotypes
replot.hapnet <- function(hapnet, labels = TRUE, scale.ratio = 20, cex = .7) {
  net <- hapnet[["net"]]
  nsamp <- hapnet[["nsamp"]]
  attr(net, "labels") <- paste0(attr(net, "labels"), " (", nsamp, ")")
  dn <- hapnet[["dn"]]
  xy <- hapnet[["xy"]]
  par_col <- c("black", "firebrick")
  par_bg <- c("grey90", "#ffd8d8")
  plot(net, xy = xy, labels = labels, scale.ratio = scale.ratio, cex = cex,
       size = nsamp, threshold = 0, show.mutation = 1,
       col = par_col[dn], bg = par_bg[dn])
  xy <- replot()
  return(xy)
}
tmp <- replot.hapnet(hapnet[["dada"]], scale.ratio = 8, cex = .7)
#Notrun hapnet[["dada"]][["xy"]] <- tmp
tmp <- replot.hapnet(hapnet[["uno5"]], scale.ratio = 8, cex = .7)
#Notrun hapnet[["uno5"]][["xy"]] <- tmp
tmp <- replot.hapnet(hapnet[["uno2"]], scale.ratio = 8, cex = .9)
#Notrun hapnet[["uno2"]][["xy"]] <- tmp
tmp <- replot.hapnet(hapnet[["mifs"]], scale.ratio = 20, cex = .4)
#Notrun hapnet[["mifs"]][["xy"]] <- tmp
rm(tmp)
# Final plotting 
hapnet.scale.ratio <- list(dada = 8, uno5 = 8, uno2 = 8, mifs = 20)
hapnet.cex <- list(data = 1.1, uno5 = 1.1, uno2 = 1.5, mifs = .6)
list(hapnet, datnam, hapnet.scale.ratio, hapnet.cex) %>%
  pmap(function(x, y, i, j) {
    svglite(paste0(path_output, "/06-Hapnet_", y, ".svg"),
            width = 16.18, height = 10)
    replot.hapnet(x, scale.ratio = i, cex = j)
    dev.off()
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
  write.csv(paste0(path_output, "/07-dn_stat.csv"))

# Save data
## Save R objects
path_saved_object <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_object, recursive = TRUE)
saveRDS(phylo_pre3, paste0(path_saved_object, "/phylo_pre3", ".obj"))
saveRDS(tree_glob, paste0(path_saved_object, "/tree_glob", ".obj"))
# Save workspace and session info
save.image(paste0(path_output, "/dnfocal.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/dnfocal.info"))
