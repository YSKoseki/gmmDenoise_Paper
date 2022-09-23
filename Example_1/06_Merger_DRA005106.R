# 06_Merger_DRA005106.R
# Last updated on 2022.4.15 by YK
# An R script to merge data from the different pipelines, DADA2, Unoise3, and MiFish
# R 4.1.2

# Packages required
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Input and output
# List of file paths
path_merger_in <- list(dada="./03-DADA2_DRA005106/05-Saved_object/seqID_tib_dada.obj",
                      uno5="./04-UNOISE3_DRA005106/G_Saved_object/seqID_tib_uno5.obj",
                      uno2="./04-UNOISE3_DRA005106/G_Saved_object/seqID_tib_uno2.obj",
                      mifs="./05-MiFishPL_DRA005106/03-Saved_object/seqID_tib_mifs.obj")
## Output directory path
path_merger_out <- "./06-Merger_DRA005106"

# Load R objects
seqID_tib_list <- lapply(path_merger_in, readRDS)

# Merge sequence ID tables from different pipelines
multi.full.join <- function(tibble.list, common.id.var="id", common.seq.var="seq") {
  library(tidyverse)
  tib.name <- names(tibble.list)
  tib.numb <- length(tibble.list)
  merged.tib <- tibble.list[[1]]
  merged.tib <- rename(merged.tib, !!tib.name[1]:=all_of(common.id.var)) # !! and := operators needed to unquote tib.name[1]
  for(i in 2:tib.numb) {
    merged.tib <- full_join(merged.tib, tibble.list[[i]], by=common.seq.var)
    merged.tib <- rename(merged.tib, !!tib.name[i]:=all_of(common.id.var))
  }
  merged.tib <- select(merged.tib, all_of(tib.name), all_of(common.seq.var))
  return(merged.tib)
}
merged_seqID_tib <- multi.full.join(seqID_tib_list)

# Give global ID to sequence 
give.glob.id <- function(sequence.id.tibble, id.var.name="glob", id.prefix="ASV_") {
  n.seqs <- nrow(sequence.id.tibble)
  seq.ids <- nchar(n.seqs) %>%
    paste0("%0", ., "d") %>% 
    paste0(id.prefix, .) %>%
    sprintf(1:n.seqs)
  new.tib <- mutate(sequence.id.tibble, !!id.var.name:=seq.ids) %>% 
    relocate(!!id.var.name)
  return(new.tib)
}
merged_seqID_tib <- give.glob.id(merged_seqID_tib)

# Make FASTA file
make.fasta2 <- function(seq.id.tibble,
                        id.var="glob",
                        seq.var="seq",
                        file.name="ASV_seqs.fa") {
  seq.id.tab <- as.data.frame(seq.id.tibble)
  id <- seq.id.tab[, id.var] %>% paste0(">", .)
  seq <- seq.id.tab[, seq.var]
  asv.seqs.out <- cbind(id, seq) %>% t() %>% c() %>% as.matrix(ncol=1)
  write.table(asv.seqs.out, file.name,
             col.names=FALSE, row.names=FALSE, quote=FALSE)
}
dir.create(path_merger_out, recursive=TRUE)
make.fasta2(seq.id.tibble=merged_seqID_tib,
            file.name=paste0(path_merger_out, "/02-ASV_seqs_merged.fa"))

# Save data
## Save R objects
path_saved_obj <- paste0(path_merger_out, "/01-Saved_object")
dir.create(path_saved_obj)
saveRDS(merged_seqID_tib, paste0(path_saved_obj, "/merged_seqID_tib", ".obj"))
## Save workspace and session info
save.image(paste0(path_merger_out, "/merger.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_merger_out, "/merger.info"))
