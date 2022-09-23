# 08_TaxonTab_DRA010703.R
# Last updated on 2022.4.16 by YK
# An R script to elaborate Claident-generated taxonomy table
# R 4.1.2

# Packages required
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Input and output
## Input R objects
path_input_seqID <- "./06-Merger_DRA010703/01-Saved_object/merged_seqID_tib.obj"
path_input_claid <- "./07-Claident_DRA010703/06-QCauto_MaxPOP05_MinSOR19+NN99.tsv"
## Output directory path
path_output <- "./08-TaxonTab_DRA010703"

# Add columns of useful taxonomic name tags
## The function
elaborate.claident.tsv <- function(path.claident.tsv, 
                                   query.id.prefix="ASV_",
                                   species.tag.var="sptag",
                                   abbrev.tag.var1="abbtag",
                                   abbrev.tag.var2="abb_id") {
  taxonomy_tib <- read_tsv(path.claident.tsv) %>%
    select(query, class, order, family, genus, species)
  # Create column of abbreviated species name
  binom <- taxonomy_tib %>% select(species) %>% unlist() %>% str_split(pattern=" ", simplify=TRUE)
  taxonomy_tib <- taxonomy_tib %>% 
    mutate(abbrev=if_else(is.na(species), NA_character_,
                          if_else(str_count(species, pattern=" ")>1, species,
                                  str_c(str_sub(binom[, 1], start=1, end=1), ". ",  binom[, 2]))))
  # Create column of species name tag, which takes any of the following forms:
  #   species (if assigned), "the-lowest-assigned-taxonomy sp.", and "UNASGD"
  taxonomy_tib <- taxonomy_tib %>% 
    mutate(!!species.tag.var:=if_else(!is.na(species), species,
                                      if_else(!is.na(genus), str_c(genus, " sp."), 
                                              if_else(!is.na(family), str_c(family, " sp."), 
                                                      if_else(!is.na(order), str_c(order, " sp."), 
                                                              if_else(!is.na(class), str_c(class, " sp."),
                                                                      "UNASGD")))))) %>% 
    # The abbreviated species name version
    mutate(!!abbrev.tag.var1:=if_else(!is.na(abbrev), abbrev,
                                      if_else(!is.na(genus), str_c(genus, " sp."), 
                                              if_else(!is.na(family), str_c(family, " sp."), 
                                                      if_else(!is.na(order), str_c(order, " sp."), 
                                                              if_else(!is.na(class), str_c(class, " sp."),
                                                                      "UNASGD")))))) %>% 
    # Create column of abbreviated species name tag with sequence ID in brackets
    mutate(!!abbrev.tag.var2:=str_c(!!sym(abbrev.tag.var1), " [", query, "]")) %>% 
    mutate(!!abbrev.tag.var2:=str_remove(!!sym(abbrev.tag.var2), pattern=query.id.prefix))
  return(taxonomy_tib)
}
## Run
claid_tib <- elaborate.claident.tsv(path_input_claid)

# Add columns of sequence IDs of different pipelines: DADA2, Unoise3, and MiFish
seqID_tib <- readRDS(path_input_seqID)
## The function
full.join.seqID <- function(seqID.tib, taxonomy.tib, by="glob", seq="seq") {
  taxonomy.tib2 <- rename(taxonomy.tib, !!by:=query)
  taxonomy.tib3 <- full_join(seqID.tib, taxonomy.tib2, by=by) %>% 
    relocate(!!seq, .after=everything())
  return(taxonomy.tib3)
}
## Run
taxon_tib <- full.join.seqID(seqID_tib, claid_tib, by="glob", seq="seq")

# Write csv file 
dir.create(path_output, recursive=TRUE)
write.csv(taxon_tib, paste0(path_output, "/02-taxon_tib.csv"))

# Save data
## Save R objects
path_saved_obj <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_obj)
saveRDS(taxon_tib, paste0(path_saved_obj, "/taxon_tib", ".obj"))
## Save workspace and session info
save.image(paste0(path_output, "/taxontab.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/taxontab.info"))
