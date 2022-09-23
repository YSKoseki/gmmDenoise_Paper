# 03_DADA2_DRA005106.R
# Last updated on 2022.4.15, by YK
# An R script to denoise amplicon sequences using the DADA2 algorithm
# For information on DADA2, see the website at https://benjjneb.github.io/dada2/index.html
# R 4.1.2

# Packages required
library(dada2); packageVersion("dada2") # 1.22.0
library(svglite); packageVersion("svglite") # 2.1.0
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Directory paths and FASTQ suffixes
## Input: FASTQ data (primer-trimmed, gzipped)
path_input_fastq <- "./02-Cutadapt_DRA005106/trimmedFASTQ"
## Input: primer-trimming process stats file
path_input_primtrim <- "./02-Cutadapt_DRA005106/primtrim.tsv"
## Output
path_dada2_out <- "./03-DADA2_DRA005106"
## Suffixes of input (primer-trimmed) and output (quality-filtered) FASTQs
fwsuff <- "_trimmed_R1.fastq.gz"
rvsuff <- "_trimmed_R2.fastq.gz"
fwsuff2 <- "_qfiltd_R1.fastq.gz"
rvsuff2 <- "_qfiltd_R2.fastq.gz"

# Plotting read quality profiles prior to quality filtering
## Get lists of input FASTQs
fwds <- list.files(path_input_fastq, pattern=fwsuff, full.names=TRUE) %>% 
  sort()
revs <- list.files(path_input_fastq, pattern=rvsuff, full.names=TRUE) %>% 
  sort()
## Get names of samples: the script below can change depending on the file name paths defined above
hiers_to_fastq <- str_count(path_input_fastq, "/") + 2
sample_names <- str_split(fwds, "/") %>% 
  sapply(`[`, hiers_to_fastq) %>% 
  str_remove(fwsuff)
## Plot
path_qprofile <- paste0(path_dada2_out, "/01-Qprofile_pre")
dir.create(path_qprofile, recursive=TRUE)
for (i in 1:length(sample_names)){
  plotQualityProfile(c(fwds[i], revs[i]))
  ggsave(paste0(sample_names[i], ".svg"), 
         path=path_qprofile, 
         device="svg", width=10, height=10, units="cm")
}
rm(i)

# Quality filtering (length filtering and abundance filtering are performed after merging the paired-end reads) 
## Preparation
path_qfiltd_fastq <- paste0(path_dada2_out, "/02-Qfiltered_FASTQ")
dir.create(path_qfiltd_fastq, recursive=TRUE)
fwds_filt <- str_split(fwds, "/") %>% 
  sapply(`[`, hiers_to_fastq) %>% 
  paste0(path_qfiltd_fastq, "/", .) %>% 
  str_replace_all(fwsuff, fwsuff2)
revs_filt <- str_replace_all(fwds_filt, "R1", "R2")
## filtering
qfiltd <- filterAndTrim(fwds, fwds_filt, 
                          revs, revs_filt,
                          maxEE=c(1, 1),
                          multithread=TRUE)
## Change the lengthy sample names in the object to the simple ones (Run IDs)
rownames(qfiltd) <- sample_names
## Plot read quality profiles
path_qprofile_qfilt <- paste0(path_dada2_out, "/03-Qprofile_post")
dir.create(path_qprofile_qfilt, recursive=TRUE)
for (i in 1:length(sample_names)){
  plotQualityProfile(c(fwds_filt[i], revs_filt[i]))
  ggsave(paste0(sample_names[i], ".svg"), 
         path=path_qprofile_qfilt, 
         device="svg", width=10, height=10, units="cm")
}
rm(i)

# Base calling error models
## The error rate learning function
set.seed(101)
fwds_error <- fwds_filt %>% 
  learnErrors(multithread=TRUE, randomize=TRUE, verbose=TRUE, 
              MAX_CONSIST=10, nbases=1e+08) # These are default values
revs_error <- revs_filt %>% 
  learnErrors(multithread=TRUE, randomize=TRUE, verbose=TRUE, 
              MAX_CONSIST=10, nbases=1e+08) # These are default values
## Base calling error profiles
path_error_plot <- paste0(path_dada2_out, "/04-Error_plot")
dir.create(path_error_plot, recursive=TRUE)
plotErrors(fwds_error, nominalQ=TRUE)
ggsave("fwds_error.svg", 
       path=path_error_plot, 
       device="svg", width=10, height=10, units="cm")
plotErrors(revs_error, nominalQ=TRUE)
ggsave("revs_error.svg", 
       path=path_error_plot, 
       device="svg", width=10, height=10, units="cm")

# Denoising
## Dereplicate sequencing reads
fwds_derep <- derepFastq(fwds_filt, verbose=TRUE)
revs_derep <- derepFastq(revs_filt, verbose=TRUE)
## Change the row names (long file names) to simple ones
names(fwds_derep) <- sample_names
names(fwds_derep) <- sample_names
## Denoise with OMEGA_A=1e-40 and OMEGA_C=1e-40 (default values)
set.seed(101)
fwds_dada <- dada(fwds_derep, err=fwds_error, 
                  OMEGA_A=1e-40, OMEGA_C=1e-40, multithread=TRUE)
revs_dada <- dada(revs_derep, err=revs_error, 
                  OMEGA_A=1e-40, OMEGA_C=1e-40, multithread=TRUE)

# Paired-end read merging
## Merge pairs of the denoised forward and reverse reads,
##   with minOverlap = 20 and maxMismatch = 0 (default values)
mergers <- mergePairs(fwds_dada, fwds_derep, 
                      revs_dada, revs_derep,
                      minOverlap=20, maxMismatch=0, verbose=TRUE)
## Make sample-by-sequence table
seq_tab <- makeSequenceTable(mergers)
getSequences(seq_tab) %>% 
  nchar() %>% 
  table()

# Additional filtering and chimera removal
## Length filtering, with the MiFish target region length ranging 163-185 bp
lfiltering <- colnames(seq_tab) %>% 
  nchar() %in% seq(163, 185)
seq_tab_lfiltd <- seq_tab[, lfiltering]
getSequences(seq_tab_lfiltd) %>% 
  nchar() %>% 
  table()
## Abundance filtering, with the cut-off value of <10 as in the MiFish Pipeline
seq_tab_nfiltd <- seq_tab_lfiltd
seq_tab_nfiltd[seq_tab_nfiltd<10] <- 0L
nonzeroseq <- colSums(seq_tab_nfiltd) %>% 
  `[`(which(.>0)) %>% 
  sort(decreasing=TRUE)
nfiltering <- colnames(seq_tab_nfiltd) %in% names(nonzeroseq)
seq_tab_nfiltd <- seq_tab_nfiltd[, nfiltering]
seq_tab_nfiltd <- seq_tab_nfiltd[, names(nonzeroseq)]
getSequences(seq_tab_nfiltd) %>% 
  nchar() %>% 
  table()
## Chimera removal  
seq_tab_nonchim <- removeBimeraDenovo(seq_tab_nfiltd, 
                                      multithread=TRUE, verbose=TRUE, method="consensus")
seq_tab_dada <- seq_tab_nonchim

# Summary table of the denoising and filtering processes
getN <- function(x) sum(getUniques(x)) 
summary_tib <- tibble(sample=sample_names,
                          ptrimd=qfiltd[, "reads.in"] %>% as.integer(),
                          qfiltd=qfiltd[, "reads.out"] %>% as.integer(),
                          dada_fwds=sapply(fwds_dada, getN) %>% as.integer(),
                          dada_revs=sapply(revs_dada, getN) %>% as.integer(),
                          merged=rowSums(seq_tab) %>% as.integer(),
                          lfiltd=rowSums(seq_tab_lfiltd) %>% as.integer(),
                          nfiltd=rowSums(seq_tab_nfiltd) %>% as.integer(),
                          nonchim=rowSums(seq_tab_nonchim) %>% as.integer())
primtrim_tib <- read_tsv(path_input_primtrim) %>% 
  transmute(sample=sample, total=as.integer(totRP))
summary_tib <- left_join(primtrim_tib, summary_tib, by="sample") %>% 
  mutate(remained=round(100*nonchim/total, 1))
summary_tib_dada <- summary_tib
write.csv(summary_tib_dada, paste0(path_dada2_out, "/06-summary_tib_dada.csv"))

# Making sequence ID table
## The function
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
## Run
seqID_tib_dada <- make.seqID.tibble(asv.table=seq_tab_dada, id.prefix="DADA")

# Making FASTA file
## The function
make.fasta <- function(asv.table,
                       header.prefix=">ASV_",
                       file.name="ASV_seqs.fa") {
  asv.seqs <- colnames(asv.table)
  n.asvs <- length(asv.seqs)
  asv.ids <- nchar(n.asvs) %>% 
    paste0("%0", ., "d") %>% 
    paste0(header.prefix, .) %>%
    sprintf(1:n.asvs)
  asv.seqs_out <- asv.ids %>% 
    rbind(asv.seqs) %>% 
    c() %>% 
    as.matrix(ncol=1)
  write.table(asv.seqs_out, file.name, 
              col.names=FALSE, row.names=FALSE, quote=FALSE)
}
## Run
make.fasta(asv.table=seq_tab_dada,
           header.prefix=">DADA",
           file.name=paste0(path_dada2_out, "/07-ASV_seqs_dada.fa"))

# Save data
## Save R objects
path_saved_obj <- paste0(path_dada2_out, "/05-Saved_object")
dir.create(path_saved_obj, recursive=TRUE)
saveRDS(seq_tab_dada, 
        paste0(path_saved_obj, "/seq_tab_dada", ".obj"))
saveRDS(summary_tib_dada, 
        paste0(path_saved_obj, "/summary_tib_dada", ".obj"))
saveRDS(seqID_tib_dada, 
        paste0(path_saved_obj, "/seqID_tib_dada", ".obj"))
## Save workspace and session info
save.image(paste0(path_dada2_out, "/dada2.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_dada2_out, "/dada2.info"))
