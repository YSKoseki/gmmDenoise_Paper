# 05_MiFishPL_DRA010703.R
# Last updated on 2022.4.15 by YK
# An R script to gather data from MiFish Pipeline outputs
# For information on MiFish Pipeline, see the website at http://mitofish.aori.u-tokyo.ac.jp/mifish/
# R 4.1.2

# Packages required
library(readxl); packageVersion("readxl") # 1.3.1
library(dada2); packageVersion("dada2") # 1.22.0
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Directory paths
## Input: MiFish Pipeline stats (xlsx) directory
path_input_stats <- "./05-MiFishPL_DRA010703/01-PL_stats"
## Input: MiFish Pipeline result (xlsx) directory
path_input_results <- "./05-MiFishPL_DRA010703/02-PL_result"
## Output
path_mifish_out <- "./05-MiFishPL_DRA010703"

# Making a summary table of MiFish Pipeline process
mifish_stats_files <- list.files(path_input_stats, pattern="xlsx", full.names=TRUE) %>% sort()
## The function
make.mifish.summary <- function(file, skip.first.rows=2) {
  mifish.summ.list <- list()
  for (i in 1:length(file)) {
    colnames <- c("sample", "fwds", "revs", "ptrimd")
    coltypes <- c("text", rep("numeric", 2), rep("skip", 14), "numeric", rep("skip",2))
    dt <- read_xlsx(file[i], sheet=1,
                    col_names=colnames,
                    col_types=coltypes,
                    skip=skip.first.rows) %>%
      mutate(remained=round(100*ptrimd/fwds, 1)) %>%  
      suppressMessages()
    mifish.summ.list[[i]] <- dt
  }
  mifish.summ <- bind_rows(mifish.summ.list) %>%
    mutate(sample=sample %>% str_remove("_"),
           fwds=fwds %>% as.integer(),
           revs=revs %>% as.integer(),
           ptrimd=ptrimd %>% as.integer(),
           remained=remained %>% as.numeric()) %>% 
    arrange(sample)
  return(mifish.summ)
}
# Run
summary_tib_mifish <- make.mifish.summary(file=mifish_stats_files)
# Write csv file
write.csv(summary_tib_mifish, paste0(path_mifish_out, "/04-summary_tib_mifs.csv"))

# Making sample-by-sequence table 
mifish_results <- list.files(path_input_results, pattern="xlsx", full.names=TRUE) %>% sort()
## The function to merge MiFish Pipeline results into one data frame
make.mifish.tibble <- function(file=NULL, skip.first.rows=0) {
  colnames <- c("sample", "species", "family", "order", "japanese", 
                "otureads", "repseq", "abundance", "confidence", "identity", 
                "lodscore", "alignlen", "mismatch", "species2", "alignlen2", 
                "mismatch2", "sequence")
  coltypes <- c("text", "text", "text", "text", "text", 
                "numeric", "text", "numeric", "text", "numeric", 
                "numeric", "numeric", "numeric", "text", "numeric", 
                "numeric", "text")
  asv.tib.list <- list()
  for (i in 1:length(file)) {
    dt <- read_xlsx(file[i],
                    sheet="BLASTN>97%Similarity",
                    col_names=colnames,
                    col_types=coltypes,
                    skip=skip.first.rows) %>% suppressWarnings()
    if (is_empty(dt)) { # Give NULL if an input file is empty
      cat(paste0("The input file [", i, "] is skipped to read because it is empty.", "\n"));
      asv.tib.list[[i]] <- NULL
    }
    else { # Filter out the header rows in the middle of data
      dt_filt <- filter(dt, sample!="Sample name")
      asv.tib.list[[i]] <- dt_filt
    }
  }
  nseg <- str_split(file, "/", simplify=TRUE) %>% ncol()
  sample.name <- str_split(file, "/", simplify=TRUE) %>% `[`(, nseg) %>% str_remove(pattern=".xlsx")
  names(asv.tib.list) <- sample.name
  asv.tib <- bind_rows(asv.tib.list)
  asv.tib$sample <- str_remove(asv.tib$sample, "_")
  return(asv.tib)
}
## Run; note the 10th file needed to handle exceptionally,
##   due to the irregular data of sample ID 'DRR069556' in the file
mifish_tib <- make.mifish.tibble(mifish_results)
## The function to make sample-by-sequence table from the 
make.sequence.table <- function(mifish.tibble, mifish.summary.tib) {
  # Split mifish.tibble by sample
  mifish.tibble <- group_by(mifish.tibble, sample)
  seq.tab.list <- select(mifish.tibble, sample, abundance, sequence) %>%
    group_split(.keep=FALSE) %>% 
    lapply(as.data.frame)
  names(seq.tab.list) <- group_keys(mifish.tibble) %>% unlist()
  # Make sample-by-sequence table with dada2::makeSequenceTable()
  mifish.seq.tab <- makeSequenceTable(seq.tab.list)
  # Add table of zero-read samples
  complete <- select(mifish.summary.tib, sample) %>% pull()
  nonzero <- rownames(mifish.seq.tab)
  zeroread <-  setdiff(complete, nonzero)
  N_samp <- length(zeroread)
  N_seq <- ncol(mifish.seq.tab)
  added.tab <- matrix(0, nrow=N_samp, ncol=N_seq, dimnames=list(zeroread, NULL))
  mifish.seq.tab2 <- rbind(mifish.seq.tab, added.tab)
  ## Sort samples
  mifish.seq.tab2 <- mifish.seq.tab2[complete, ]
  return(mifish.seq.tab2)
}
## Run
seq_tab_mifish <- make.sequence.table(mifish.tibble=mifish_tib,
                                      mifish.summary.tib=summary_tib_mifish)

# Making sequence ID table
## The function
make.seqID.tibble <- function(asv.table, id.prefix="ASV") {
  asv.seqs <- colnames(asv.table)
  n.asvs <- length(asv.seqs)
  asv.ids <- nchar(n.asvs) %>% 
    paste0("_", "%0", ., "d") %>% 
    paste0(id.prefix, .) %>%
    sprintf(1:n.asvs)
  seqID.tib <- tibble(id=as.character(asv.ids), seq=as.character(asv.seqs))
  return(seqID.tib)
}
## Run
seqID_tib_mifish <- make.seqID.tibble(seq_tab_mifish, id.prefix="MIFS")

# Making FASTA file
make.fasta <- function(asv.table,
                       header.prefix=">ASV_",
                       file.name="ASV_seqs.fa") {
  asv.seqs <- colnames(asv.table)
  n.asvs <- length(asv.seqs)
  asv.ids <- nchar(n.asvs) %>% 
    paste0("%0", ., "d") %>% 
    paste0(header.prefix, .) %>%
    sprintf(1:n.asvs)
  asv.seqs.out <- asv.ids %>% 
    rbind(asv.seqs) %>% 
    c() %>% 
    as.matrix(ncol=1)
  write.table(asv.seqs.out, file.name, 
              col.names=FALSE, row.names=FALSE, quote=FALSE)
}
make.fasta(asv.table=seq_tab_mifish,
           header.prefix=">MIFS",
           file.name=paste0(path_mifish_out, "/05-ASV_seqs_mifs.fa"))

# Save data
## Save R objects
path_saved_obj <- paste0(path_mifish_out, "/03-Saved_object")
dir.create(path_saved_obj)
saveRDS(summary_tib_mifish, paste0(path_saved_obj, "/summary_tib_mifs", ".obj"))
saveRDS(seq_tab_mifish, paste0(path_saved_obj, "/seq_tab_mifs", ".obj"))
saveRDS(seqID_tib_mifish, paste0(path_saved_obj, "/seqID_tib_mifs", ".obj"))
## Save workspace and session info
save.image(paste0(path_mifish_out, "/mifishpl.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_mifish_out, "/mifishpl.info"))
