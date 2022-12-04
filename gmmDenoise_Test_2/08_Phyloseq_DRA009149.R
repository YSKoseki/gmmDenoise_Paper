# 08_Phyloseq_DRA009149.R
# Last updated on 2022.11.13 by YK
# An R script to construct phyloseq objects with data from different pipelines
# R 4.1.2

# Packages required
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9018
library(Biostrings); packageVersion("Biostrings") # 2.62.0
library(DECIPHER); packageVersion("DECIPHER") # 2.22.0
library(phangorn); packageVersion("phangorn") # 2.8.1
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Input and output
## Input R objects: 
##   sequence-read-abundance tables (from different pipelines), taxonomy table, and sample table
path_seq_tab <- list(dada="./03-DADA2_DRA009149/05-Saved_object/seq_tab_dada.obj",
                     uno3="./04-UNOISE3_DRA009149/K_Saved_object/seq_tab_uno3.obj",
                     nodn="./04-UNOISE3_DRA009149/K_Saved_object/seq_tab_nodn.obj")
path_taxon_tib <- "./07-TaxonTab_DRA009149/01-Saved_object/taxon_tib.obj"
path_samp_tab <- "./08-SraRunTable.csv"
## Output directory
path_output <- "./08-Phyloseq_DRA009149"

# Construct a series of phyloseq objects (stored in one list)
## The function
construct.phyloseq <- function(path.seq.tab, path.taxon.tib, path.samp.tab,
                               taxa.are.rows=FALSE,
                               seq.var="seq", samp.var="Run") {
  # Sequence read abundance data
  seq.tab <- readRDS(path.seq.tab) %>% otu_table(taxa_are_rows=taxa.are.rows)
  # Taxonomy data
  taxon.tab <- readRDS(path.taxon.tib) %>% column_to_rownames(seq.var) %>%
    as.matrix() %>% tax_table()
  # Sample data
  samp.tab <- read.csv(path_samp_tab) %>% column_to_rownames(samp.var) %>%
    sample_data()
  # Construct phyloseq object
  ps <- phyloseq(seq.tab, taxon.tab, samp.tab)
  return(ps)
}
## Run
phylo <- pmap(
  list(path_seq_tab, path_taxon_tib, path_samp_tab),
  function(x, y, z) {
    construct.phyloseq(path.seq.tab=x, path.taxon.tib=y, path.samp.tab=z)
  }
)

# Use ASV IDs (the global) instead of full DNA sequences as sequence identifiers for convenience
## The function
fullseq.to.id <- function(physeq, id.var="abb_id") {
  dna <- taxa_names(physeq) %>% Biostrings::DNAStringSet()
  names(dna) <- taxa_names(physeq)
  physeq <- merge_phyloseq(physeq, dna)
  taxa_names(physeq) <- tax_table(physeq) %>% `[`(, id.var) %>% as.vector()
  return(physeq)
}
## Run
phylo <- lapply(phylo, fullseq.to.id); phylo

# Build phylogenetic trees and bind them to the phyloseq object
## Multiple alignment with DECIPHER
##   For 'DECIPHER', see the document at https://rdrr.io/bioc/DECIPHER/f/inst/doc/ArtOfAlignmentInR.pdf
alg <- phylo %>% 
  lapply(refseq) %>% 
  lapply(function(x) DECIPHER::AlignSeqs(x, processors=NULL)) %>% 
  lapply(function(x) DECIPHER::StaggerAlignment(x, processors=NULL))
## Compute pairwise distances of aligned sequences
##   For 'phangorn', see vignette("Trees", package="phangorn")
phydat <- alg %>% 
  lapply(as.matrix) %>% 
  lapply(phangorn::phyDat)
phydat_dist <- phydat %>% 
  lapply(phangorn::dist.ml)
## Build neighbor-Joining tree
treeNJ <- phydat_dist %>% 
  lapply(phangorn::NJ)
##Build Jukes-Cantor model ML tree
treeJC <- map2(treeNJ, phydat, function(x, y) phangorn::pml(x, data=y)) %>%
  lapply(function(x) optim.pml(x, optNni=TRUE))
## Build GTR+I+G model ML tree (General Time-Reversible model, with corrections for Invariant characters and Gamma-distributed rate heterogeneity)
treeGTR <- map2(treeNJ, phydat, function(x, y) phangorn::pml(x, data=y)) %>%
  lapply(function(x) update(x, k=4, inv=0.2)) %>%
  lapply(function(x) phangorn::optim.pml(x, model="GTR", 
                                         optInv=TRUE, optGamma=TRUE,
                                         rearrangement="stochastic",
                                         control=phangorn::pml.control(trace=0)
                                         ))
## Model selection with AIC
mod.sel <- rbind(treeJC=lapply(treeJC, AIC) %>% unlist(),
                 treeGTR=lapply(treeGTR, AIC) %>% unlist()); mod.sel
## Bind the tree to the phyloseq object
phylo <- map2(phylo, treeGTR,
              function(x, y) merge_phyloseq(x, y$tree))

# Save data
## Save R objects
path_saved_object <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_object, recursive=TRUE)
saveRDS(phylo, paste0(path_saved_object, "/phylo", ".obj"))
## Save workspace and session info
save.image(paste0(path_output, "/phyloseq.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/phyloseq.info"))
