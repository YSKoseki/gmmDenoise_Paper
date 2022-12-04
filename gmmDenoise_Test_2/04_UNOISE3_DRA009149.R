# 04_UNOISE3_DRA009149.R
# Last updated on 2022.11.11 by YK
# An R script to extract amplicon sequence variants with and without denoising
#   using the UNOISE3 algorithms in the JAMP package
# For information on JAMP, see the website at https://github.com/VascoElbrecht/JAMP
# R 4.1.2

# Packages required
library(R.utils); packageVersion("R.utils") # 2.11.0
library(JAMP); packageVersion("JAMP") # 0.77
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Directory paths and FASTQ suffixes
## Input FASTQ data (primer-trimmed, gzipped)
path_input_fastq <- "./02-Cutadapt_DRA009149/trimmedFASTQ"
## Output
path_UNOISE3_out <- "./04-UNOISE3_DRA009149"
## Suffixes of input (primer-trimmed) FASTQs
fwsuff <- "_trimmed_R1.fastq.gz"
rvsuff <- "_trimmed_R2.fastq.gz"

# A. Preprocess for the JAMP analysis
## Copy and unzip the FASTQ gz files
fwds <- list.files(path_input_fastq, pattern=fwsuff, full.names=TRUE) %>% 
  sort()
revs <- list.files(path_input_fastq, pattern=rvsuff, full.names=TRUE) %>% 
  sort()
dir.create("_data", recursive=TRUE)
file.copy(from=c(fwds, revs), to="_data", recursive=TRUE)
list.files("_data", full.names=TRUE) %>% 
  for (i in .) {gunzip(i)}
## Create JAMP working directory and move the unzipped FASTQs there 
path_prev_wd <- getwd()
dir.create(path_UNOISE3_out, recursive=TRUE)
setwd(path_UNOISE3_out)
path_curr_wd <- getwd() 
Empty_folder() # Creates "./A_Empty_Folder/_data and _stats")
file.copy(from=paste0(path_prev_wd, "/_data"), to="./A_Empty_Folder/", recursive=TRUE)
## Remove the source files, then the source directory
list.files(paste0(path_prev_wd, "/_data"), full.names=TRUE) %>%
  file.remove() 
file.remove(paste0(path_prev_wd, "/_data"))

# B. Merging paired-end reads with fastq_maxdiffs=0, fastq_maxdiffpct=0 and fastq_minovlen=12
#   (so as to match with DADA2) 
Merge_PE(fastq_maxdiffs=0, fastq_maxdiffpct=0, fastq_minovlen=12)

# C. Length filtering, set min length at 100 bp (Tsuji et al. 2020)
Minmax(min=100, max=1000)

# D. Quality filtering
Max_ee(max_ee=1)

# E. Denoising with alpha=3 and
#   abundance filtering with minsize=4 (following Tsuji et al. 2020)
list.files("D_U_max_ee/_data", full.names=TRUE) %>% 
  Denoise(strategy="unoise3", unoise_alpha=3, minsize=4, 
          minrelsize=0, OTUmin=0, minhaplosize=0, withinOTU=0) # No additional filtering

# F. A virtual no denoising, with alpha=10000 and
#   abundance filtering with minsize=4 (following Tsuji et al. 2020)
list.files("D_U_max_ee/_data", full.names=TRUE) %>% 
  Denoise(strategy="unoise3", unoise_alpha=10000, minsize=4, 
          minrelsize=0, OTUmin=0, minhaplosize=0, withinOTU=0) # No additional filtering

# G. Data assembling for downstream analyses
## Making a summary table of the denoising and filtering processes
### The function
make.summary.tibble <- function(ptrim.stats.tsv=NULL,
                                merge.stats.csv=NULL,
                                lfiltd.stats.csv=NULL,
                                qfiltd.stats.csv=NULL,
                                unoise.stats.txt=NULL,
                                input.fastq.suff="_[rR][12].fastq.gz") {
  namesuff<- str_remove(input.fastq.suff, "_[rR][12].fastq.gz")
  namesuff2 <- str_replace(input.fastq.suff, "_[rR][12].fastq.gz", "_PE_minmax.fastq")
  # Table of primer-trimming process
  ptrim.stats <- read_tsv(ptrim.stats.tsv) %>% 
    transmute(sample=sample,
              total=totRP %>% as.integer(),
              ptrimd=filtdRP %>% as.integer()) %>% 
    suppressMessages()
  # Table of paired-end merging process
  merge.stats <- read_csv(merge.stats.csv) %>% 
    transmute(sample=str_remove_all(Sample, namesuff),
              ptrimd=Sequ_count_in %>% as.integer(), 
              merged=Sequ_count_out %>% as.integer()) %>% 
    suppressMessages()
  # Table of length-filtering process
  lfiltd.stats <- read_csv(lfiltd.stats.csv) %>%
    transmute(sample=str_remove_all(Sample, namesuff2),
              merged=Sequ_count_in %>% as.integer(),
              lfiltd=Sequ_count_out %>% as.integer()) %>%
    suppressMessages()
  # Table of quality-filtering process
  qfiltd.stats <- read_csv(qfiltd.stats.csv) %>%
    transmute(sample=str_remove_all(Sample, namesuff),
              lfiltd=Sequ_count_in %>% as.integer(),
              qfiltd=Sequ_count_out %>% as.integer()) %>% 
    suppressMessages()
  # Table of UNOISE3-denoising process
  unoise.stats.txt <- file(unoise.stats.txt, open="r")
  samplename.list <- list()
  unoise.list <- list()
  nonchim.list <- list()
  i <- 0
  j <- 0
  while (TRUE) {
    line <- readLines(unoise.stats.txt, n=1)
    if (length(line)==0) {
      break
    }
    # Sample names
    namesuff<- namesuff
    targetline1 <- str_detect(line, "Reading file") &&
      str_detect(line, "fasta")
    if (targetline1) {
      i <- i+1
      sample <- str_split(line, "1_derep/", simplify=TRUE) %>% 
        `[`(2) %>%
        str_split(namesuff, simplify=TRUE) %>% 
        `[`(1)
      samplename.list[[i]] <- sample
    }
    # Stats
    targetline2 <- str_detect(line, "Found", negate=TRUE) &&
      str_detect(line, "chimeras") &&
      str_detect(line, "non-chimeras")
    if (targetline2) {
      j <- j+1
      seg <- str_split(line, ", ", simplify=TRUE)
      chimera <- parse_number(seg[1])
      nonchim <- parse_number(seg[2])
      unoised <- chimera + nonchim
      unoise.list[[j]] <- unoised
      nonchim.list[[j]] <- nonchim
    }
  }
  unoise.stats <- tibble(sample=unlist(samplename.list), 
                         unoise=unlist(unoise.list) %>% as.integer(),
                         nonchim=unlist(nonchim.list) %>% as.integer())
  summary.tib <-  ptrim.stats %>% 
    left_join(select(merge.stats, -ptrimd), by="sample") %>% 
    left_join(select(lfiltd.stats, -merged), by="sample") %>% 
    left_join(select(qfiltd.stats, -lfiltd), by="sample") %>% 
    left_join(unoise.stats, by="sample") %>% 
    mutate(remained=round(nonchim/total*100, 1)) %>% 
    replace_na(list(total=0, ptrimd=0, merged=0, lfiltd=0, qfiltd=0, unoise=0, nonchim=0, remained=0))
  close(unoise.stats.txt)
  return(summary.tib)
}
### Tables of the analyses
summary_tib_uno3 <- make.summary.tibble(ptrim.stats.tsv="../02-Cutadapt_DRA009149/primtrim.tsv",
                                        merge.stats.csv="./B_merge_PE/_stats/B_sequ_length_abund.csv",
                                        lfiltd.stats.csv="./C_Minmax/_stats/C_minmax_pass.csv",
                                        qfiltd.stats.csv="./D_U_max_ee/_stats/D_max_ee_stats.csv",
                                        unoise.stats.txt="./E_Denoising/_stats/2_denoise_logs.txt",
                                        input.fastq.suff="_trimmed_R1.fastq.gz")
write.csv(summary_tib_uno3, "./G_summary_tib_uno3.csv")

summary_tib_nodn <- make.summary.tibble(ptrim.stats.tsv="../02-Cutadapt_DRA009149/primtrim.tsv",
                                        merge.stats.csv="./B_merge_PE/_stats/B_sequ_length_abund.csv",
                                        lfiltd.stats.csv="./C_Minmax/_stats/C_minmax_pass.csv",
                                        qfiltd.stats.csv="./D_U_max_ee/_stats/D_max_ee_stats.csv",
                                        unoise.stats.txt="./F_Denoising/_stats/2_denoise_logs.txt",
                                        input.fastq.suff="_trimmed_R1.fastq.gz")
write.csv(summary_tib_nodn, "./H_summary_tib_nodn.csv")

## Making ASV table (sample-by-sequence matrix of amplicon reads)
### The function
make.seq.table <- function(haplo.table.file) {
  asv_tib <- read_csv(haplo.table.file) %>% 
    filter(!is.na(haplotype)) %>% 
    suppressMessages()
  haplo_id <- select(asv_tib, haplotype) %>% as.matrix() %>% parse_number()
  asv_id <- nrow(asv_tib) %>%
    nchar(.) %>%
    paste0("%0", ., "d") %>%
    sprintf(haplo_id) %>% 
    paste0("asv_", .)
  asv_tib <- mutate(asv_tib, asv=asv_id) %>% 
    select(-sort, -haplotype, -OTU) %>% 
    select(asv, sequences, everything()) %>% 
    arrange(asv)
  colnames(asv_tib) %<>% 
    str_split("_", simplify=TRUE) %>% 
    `[`(, 1)
  asv_tib <- select(asv_tib, -asv) %>% 
    column_to_rownames("sequences") %>% 
    t()
}
### Run
seq_tab_uno3 <- make.seq.table("./E_Denoising/E_haplo_table.csv")
seq_tab_nodn <- make.seq.table("./F_Denoising/E_haplo_table.csv")

## Making sequence ID table
### The function
make.seqID.tibble <- function(asv.table=NULL,
                              id.prefix="ASV") {
  asv.seqs <- colnames(asv.table)
  n.asvs <- length(asv.seqs)
  asv.ids <- nchar(n.asvs) %>% 
    paste0("_", "%0", ., "d") %>% 
    paste0(id.prefix, .) %>%
    sprintf(1:n.asvs)
  seqID.tib <- tibble(id=as.character(asv.ids), seq=as.character(asv.seqs))
  return(seqID.tib)
}
### Run
seqID_tib_uno3 <- make.seqID.tibble(asv.table=seq_tab_uno3, id.prefix="UNO3")
seqID_tib_nodn <- make.seqID.tibble(asv.table=seq_tab_nodn, id.prefix="NODN")

## Making FASTA file
### The function
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
### Run
make.fasta(asv.table=seq_tab_uno3,
           header.prefix=">UNO3",
           file.name="I_ASV_seqs_uno3.fa")
make.fasta(asv.table=seq_tab_nodn,
           header.prefix=">NODN",
           file.name="J_ASV_seqs_nodn.fa")

## Save data
## Save R objects
setwd(path_curr_wd)
path_saved_obj <- "./K_Saved_object"
dir.create(path_saved_obj, recursive=TRUE)
saveRDS(seq_tab_uno3, 
        paste0(path_saved_obj, "/seq_tab_uno3", ".obj"))
saveRDS(seq_tab_nodn, 
        paste0(path_saved_obj, "/seq_tab_nodn", ".obj"))
saveRDS(summary_tib_uno3, 
        paste0(path_saved_obj, "/summary_tib_uno3", ".obj"))
saveRDS(summary_tib_nodn, 
        paste0(path_saved_obj, "/summary_tib_nodn", ".obj"))
saveRDS(seqID_tib_uno3, 
        paste0(path_saved_obj, "/seqID_tib_uno3", ".obj"))
saveRDS(seqID_tib_nodn, 
        paste0(path_saved_obj, "/seqID_tib_nodn", ".obj"))
### Save workspace and session info
save.image(paste0("unoise3.RData"))
writeLines(capture.output(sessionInfo()), paste0("unoise3.info"))
