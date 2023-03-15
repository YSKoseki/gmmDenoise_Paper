# 12_DNsmry_DRA010703.R
# Last updated on 2023.3.15 by YK
# An R script to evaluate denoising effects
# R 4.1.2

# Packages required
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9018
library(cowplot); packageVersion("cowplot") # 1.1.1
library(ggtree); packageVersion("ggtree") # 3.2.1
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Paths
## List of phyloseq objects
path_input1 <- "./10-Phyloseq_DRA010703/01-Saved_object/phylo.obj"
## GMM-inferred ASV read count cutoff threshold
path_input2 <- "./11-GMMdn_DRA010703/01-Saved_object/thresh_tab.obj"
## Output directory
path_output <- "./12-DNsmry_DRA010703"

# Load the phyloseq objects
phylo_pre <- readRDS(path_input1)
thresh_tab <- readRDS(path_input2)

# Create denoised ASV data sets
thresh_lis <- thresh_tab %>%
  select(data, norm) %>%
  pivot_wider(names_from = data, values_from = norm) %>% 
  as.list() 
phylo_post <- thresh_lis %>% 
  map2(phylo_pre, function(x, y) {
    # In case of cutoff value not available (no denoising)
    if (is.na(x)) y
    # Retain ASVs with a read count larger than the cutoff
    else prune_taxa(taxa_sums(y) > x, y)
  })

# Summary stats for denoising effects
## Function to calculate summary stats
smry.stats <- function(phydat.pre, phydat.post) {
  # Number of ASVs
  asv <- tibble(pre = phyloseq::ntaxa(phydat.pre), 
                post = phyloseq::ntaxa(phydat.post),
                retain = post / pre)
  # Number of ASVs by class
  pre_tib <- phydat.pre %>%
    tax_table() %>% as.data.frame() %>% rownames_to_column(var = "id")
  post_tib <- phydat.post %>%
    tax_table() %>% as.data.frame() %>% rownames_to_column(var = "id")
  asv_cl_pre <- pre_tib %>% group_by(class) %>% summarize(pre = n())
  asv_cl_post <- post_tib %>% group_by(class) %>% summarize(post = n())
  asv_by_cl <- full_join(asv_cl_pre, asv_cl_post, by = "class") %>%
    arrange(desc(post)) %>% mutate(retain = post / pre)
  # Number of identified species
  sp_pre <- pre_tib %>%
    select(class, species) %>% filter(species != "NA_character_") %>%
    distinct(species, .keep_all = TRUE)
  sp_post <- post_tib %>%
    select(class, species) %>% filter(species != "NA_character_") %>%
    distinct(species, .keep_all = TRUE)
  sp <- tibble(pre = sp_pre %>% nrow(),
               post = sp_post %>% nrow(),
               retain = post / pre)
  # Number of identified species by class
  sp_cl_pre <- sp_pre %>% group_by(class) %>% summarize(pre = n())
  sp_cl_post <- sp_post %>% group_by(class) %>% summarize(post = n())
  sp_by_cl <- full_join(sp_cl_pre, sp_cl_post, by = "class") %>%
    arrange(desc(post)) %>% mutate(retain = post / pre)
  # Fish species
  fsclass <- c("Actinopteri", "Chondrichthyes")
  fs_pre <- sp_pre %>% filter(class %in% fsclass) %>%
    select(species) %>% group_by(species) %>% summarize(pre = n())
  fs_post <- sp_post %>% filter(class %in% fsclass) %>%
    select(species) %>% group_by(species) %>% summarize(post = n())
  fs_sp <- full_join(fs_pre, fs_post, by = "species") %>%
    arrange(desc(post))
  # Non-fish species
  nonfs_pre <- sp_pre %>% filter(!(class %in% fsclass)) %>%
    select(species) %>% group_by(species) %>% summarize(pre = n())
  nonfs_post <- sp_post %>% filter(!(class %in% fsclass)) %>%
    select(species) %>% group_by(species) %>% summarize(post = n())
  nonfs_sp <- full_join(nonfs_pre, nonfs_post, by = "species") %>%
    arrange(desc(post))
  # List object
  l <- list(asv = asv, asv_by_cl = asv_by_cl,
            sp = sp, sp_by_cl = sp_by_cl,
            fs_sp = fs_sp, nonfs_sp = nonfs_sp)
  return(l)
}
## Run
smrystats <- list(phylo_pre, phylo_post) %>%
  pmap(smry.stats)

# Read size by ASV before and after denoising
## Function
asv.lists <- function(phydat.pre, phydat.post) {
  # Taxonomy and abundance data
  melt_pre <- psmelt(phydat.pre) %>%
    select(OTU, class, order, family, genus, species, glob, Abundance) %>%
    group_by(OTU) %>%
    summarize(class = first(class), order = first(order),
              family = first(family), genus = first(genus),
              species = first(species), asv = first(glob),
              pre = sum(Abundance))
  melt_post <- psmelt(phydat.post) %>%
    select(OTU, class, order, family, genus, species, glob, Abundance) %>%
    group_by(OTU) %>%
    summarize(class = first(class), order = first(order),
              family = first(family), genus = first(genus),
              species = first(species), asv = first(glob),
              post = sum(Abundance)) 
  # Sequence data
  seq_pre <- refseq(phydat.pre) %>% as.data.frame() %>%
    rownames_to_column() %>% rename("OTU" = rowname, "sequence" = x)
  seq_post <- refseq(phydat.post) %>% as.data.frame() %>%
    rownames_to_column() %>% rename("OTU" = rowname, "sequence" = x)
  # Combined data
  join_pre <- full_join(melt_pre, seq_pre, by = "OTU") %>% select(-OTU)
  join_post <- full_join(melt_post, seq_post, by = "OTU") %>% select(-OTU)
  join_dat <- full_join(join_pre, join_post,
                        by = c("class", "order", "family", "genus",
                               "species", "asv", "sequence")) %>%
    select(class, order, family, genus, species, asv, sequence, pre, post) %>%
    arrange(class, order, family, genus, species, asv)
  # Fish data
  fsclass <- c("Actinopteri", "Chondrichthyes")
  fs_dat <- join_dat %>% filter(class %in% fsclass)
  # Non-fish data
  nonfs_dat <- join_dat %>% filter(!(class %in% fsclass))
  # List object
  l <- list(fs_asv = fs_dat, nonfs_asv = nonfs_dat)
  return(l)
}
## Run
asvlist <- list(phylo_pre, phylo_post) %>%
  pmap(asv.lists)

## Save stats
dir.create(path_output, recursive=TRUE)
### Number of ASV
smrystats %>% 
  lapply(function(x) x[["asv"]]) %>% 
  bind_rows() %>%
  mutate(data = names(smrystats)) %>% 
  select(data, everything()) %>%
  write.csv(paste0(path_output, "/02-diff_asv.csv"))
### Number of ASVs by class
smrystats %>%
  lapply(function(x) x[["asv_by_cl"]]) %>% 
  map2(list("dada", "uno5", "uno2", "mifs"),
       function(x, y) {
         mutate(x, data = y) %>% select(data, everything())
       }) %>%
  bind_rows() %>%
  write.csv(paste0(path_output, "/03-diff_asv_by_cl.csv"))
### Number of identified species
smrystats %>% 
  lapply(function(x) x[["sp"]]) %>% 
  bind_rows() %>%
  mutate(data = names(smrystats)) %>% 
  select(data, everything()) %>%
  write.csv(paste0(path_output, "/04-diff_sp.csv"))
### Number of identified species by class
smrystats %>% 
  lapply(function(x) x[["sp_by_cl"]]) %>% 
  map2(list("dada", "uno5", "uno2", "mifs"),
       function(x, y) {
         mutate(x, data = y) %>% select(data, everything())
       }) %>%
  bind_rows() %>%
  write.csv(paste0(path_output, "/05-diff_sp_by_cl.csv"))
## Identification of fish species
temp <- smrystats %>% 
  lapply(function(x) x[["fs_sp"]])
temp[["dada"]] %>% 
  full_join(temp[["uno5"]], suffix = c("", ".uno5"), by = "species") %>%
  full_join(temp[["uno2"]], suffix = c("", ".uno2"), by = "species") %>%
  full_join(temp[["mifs"]], suffix = c("", ".mifs"), by = "species") %>%
  mutate(pre.dada = pre, post.dada = post, .keep = "unused") %>%
  select(species, pre.dada, post.dada, everything()) %>%
  write.csv(paste0(path_output, "/06-diff_fs_sp.csv"))
rm(temp)
## Identification of non-fish species
temp2 <- smrystats %>% 
  lapply(function(x) x[["nonfs_sp"]])
temp2[["dada"]] %>% 
  full_join(temp2[["uno5"]], suffix = c("", ".uno5"), by = "species") %>%
  full_join(temp2[["uno2"]], suffix = c("", ".uno2"), by = "species") %>%
  full_join(temp2[["mifs"]], suffix = c("", ".mifs"), by = "species") %>%
  mutate(pre.dada = pre, post.dada = post, .keep = "unused") %>%
  select(species, pre.dada, post.dada, everything()) %>%
  write.csv(paste0(path_output, "/07-diff_nonfs_sp.csv"))
rm(temp2)
## List of all detected fish ASVs
temp3 <- asvlist %>%
  lapply(function(x) x[["fs_asv"]])
temp3[["dada"]] %>% 
  full_join(temp3[["uno5"]], suffix = c("", ".uno5"),
            by = c("class", "order", "family", "genus",
                   "species", "asv", "sequence")) %>%
  full_join(temp3[["uno2"]], suffix = c("", ".uno2"),
            by = c("class", "order", "family", "genus",
                   "species", "asv", "sequence")) %>%
  full_join(temp3[["mifs"]], suffix = c("", ".mifs"),
            by = c("class", "order", "family", "genus",
                   "species", "asv", "sequence")) %>%
  mutate(pre.dada = pre, post.dada = post, .keep = "unused") %>%
  select(class, order, family, genus, species, asv, sequence,
         pre.dada, post.dada, everything()) %>%
  arrange(class, order, family, genus, species, asv) %>% 
  write.csv(paste0(path_output, "/08-diff_fs_asv.csv"))
rm(temp3)
## List of all detected non-fish ASVs
temp4 <- asvlist %>%
  lapply(function(x) x[["nonfs_asv"]])
temp4[["dada"]] %>% 
  full_join(temp4[["uno5"]], suffix = c("", ".uno5"),
            by = c("class", "order", "family", "genus",
                   "species", "asv", "sequence")) %>%
  full_join(temp4[["uno2"]], suffix = c("", ".uno2"),
            by = c("class", "order", "family", "genus",
                   "species", "asv", "sequence")) %>%
  full_join(temp4[["mifs"]], suffix = c("", ".mifs"),
            by = c("class", "order", "family", "genus",
                   "species", "asv", "sequence")) %>%
  mutate(pre.dada = pre, post.dada = post, .keep = "unused") %>%
  select(class, order, family, genus, species, asv, sequence,
         pre.dada, post.dada, everything()) %>%
  arrange(class, order, family, genus, species, asv) %>% 
  write.csv(paste0(path_output, "/09-diff_nonfs_asv.csv"))
rm(temp4)

# Save data
## Save R objects
path_saved_object <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_object, recursive=TRUE)
saveRDS(phylo_pre, paste0(path_saved_object, "/phylo_pre", ".obj"))
saveRDS(phylo_post, paste0(path_saved_object, "/phylo_post", ".obj"))
# Save workspace and session info
save.image(paste0(path_output, "/dnsmry.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/dnsmry.info"))
