# 02_sequencing.R
# A simulation of DNA sequencing of PCR products

# Packages required
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Parameters 
repl_n <- 100
cyc <- 35
tot_read_n <- 1.5e+06

# Input directories and file name patterns
path_in <- "./01-pcr.R"
path_repl <- paste0(path_in, "/repl_",
                    formatC(1:repl_n, width = nchar(repl_n), flag = "0"))
patrn_prod <- formatC(0:cyc, width = nchar(cyc), flag = "0") %>%
  paste0("products_", .) 
patrn_read <- "read_size"

# Output directories
path_out <- "./02-sequencing.R"
dir.create(path_out, recursive = TRUE)
path_obj <- paste0(path_out, "/01-Saved_obj")
dir.create(path_obj, recursive = TRUE)

# Read PCR product data into one data frame (tibble)
products <- path_repl %>%
  list.files(pattern = patrn_prod[cyc + 1], full.names = TRUE) %>% sort() %>%
  lapply(function(x) {
    read_csv(x, col_select = c(repl, cyc, seq, conc), show_col_types = FALSE)
  }) %>%
  lapply(function(x) {
    y <- x %>%
      mutate(
        seq_type = if_else(str_detect(seq, pattern = "e$"), "False", "True")
      ) %>%
      select(repl, cyc, seq_type, seq, conc)
    return(y)
  }) %>%
  bind_rows() %>% group_nest(repl)

# Read "read size" data into one data frame (tibble)
read_size <- path_repl %>%
  list.files(pattern = patrn_read, full.names = TRUE) %>% sort() %>%
  lapply(function(x) {
    y <- x %>%
      read_csv(
        col_select = c(repl, cyc, N, sum_conc),
        show_col_types = FALSE
      ) %>%
      filter(cyc == 35)
    return(y)
  }) %>%
  bind_rows()

# Simulation of sequencing: a multinomial sampling of sequencing reads
set.seed(0829)
sequencing <- map2(
  as.list(products$data),
  as.list(read_size$sum_conc),
  function(x, y) {
    z <- x %>%
      mutate(
        reads = rmultinom(1, size = tot_read_n, prob = conc / y) %>% c()
      ) %>%
      filter(reads > 0)
    return(z)
  }
) %>%
  bind_rows(.id = "repl") %>%
  mutate(repl = as.numeric(repl))

# Summary of simulated sequencing
## By sequence type (true and error sequences)
summary_tib <- sequencing %>% nest_by(repl, seq_type) %>%
  mutate(
    seq_n = length(data$seq),
    sum_reads = sum(data$reads)
  ) %>%
  select(repl, seq_type, seq_n, sum_reads)
write_csv(summary_tib, paste0(path_out, "/02-summary_tib.csv"))

summary_seq_tib <- summary_tib %>% nest_by(seq_type) %>%
  mutate(
    mean_seq_n = mean(data$seq_n),
    sd_seq_n = sd(data$seq_n),
    min_seq_n = min(data$seq_n),
    max_seq_n = max(data$seq_n)
  ) %>%
  select(seq_type, mean_seq_n, sd_seq_n, min_seq_n, max_seq_n)
write_csv(summary_seq_tib, paste0(path_out, "/03-summary_seq_tib.csv"))

summary_sumread_tib <- summary_tib %>% nest_by(seq_type) %>%
  mutate(
    mean_sum_reads = mean(data$sum_reads),
    sd_sum_reads = sd(data$sum_reads),
    min_sum_reads = min(data$sum_reads),
    max_sum_reads = max(data$sum_reads)
  ) %>%
  select(seq_type, mean_sum_reads, sd_sum_reads, min_sum_reads, max_sum_reads)
write_csv(summary_sumread_tib, paste0(path_out, "/04-summary_sumread_tib.csv"))

summary_seqread_tib <- sequencing %>% nest_by(seq_type) %>%
  mutate(
    mean_read_n = mean(data$reads),
    sd_read_n = sd(data$reads),
    min_read_n = min(data$reads),
    max_read_n = max(data$reads)
  ) %>%
  select(seq_type, mean_read_n, sd_read_n, min_read_n, max_read_n)
write_csv(summary_seqread_tib, paste0(path_out, "/05-summary_seqread_tib.csv"))

## By error type (single- and double-base-substituted sequences)
summary_err_tib <- sequencing %>%
  filter(seq_type == "False") %>%
  mutate(nsub = seq %>% str_count("e") %>% `-`(1) %>% factor()) %>%
  nest_by(repl, nsub, .drop = FALSE) %>%
  mutate(
    seq_n = length(data$seq),
    sum_reads = sum(data$reads)
  ) %>%
  select(repl, nsub, seq_n, sum_reads)
write_csv(summary_err_tib, paste0(path_out, "/06-summary_err_tib.csv"))

summary_err_seq_tib <- summary_err_tib %>% nest_by(nsub) %>%
  mutate(
    mean_seq_n = mean(data$seq_n),
    sd_seq_n = sd(data$seq_n),
    min_seq_n = min(data$seq_n),
    max_seq_n = max(data$seq_n)
  ) %>%
  select(nsub, mean_seq_n, sd_seq_n, min_seq_n, max_seq_n)
write_csv(summary_err_seq_tib, paste0(path_out, "/07-summary_err_seq_tib.csv"))

summary_err_sumread_tib <- summary_err_tib %>% nest_by(nsub) %>%
  mutate(
    mean_sum_reads = mean(data$sum_reads),
    sd_sum_reads = sd(data$sum_reads),
    min_sum_reads = min(data$sum_reads),
    max_sum_reads = max(data$sum_reads)
  ) %>%
  select(nsub, mean_sum_reads, sd_sum_reads, min_sum_reads, max_sum_reads)
write_csv(summary_err_sumread_tib, paste0(path_out, "/08-summary_err_sumread_tib.csv"))

summary_err_seqread_tib <- sequencing %>%
  filter(seq_type == "False") %>%
  mutate(nsub = seq %>% str_count("e") %>% `-`(1) %>% factor()) %>%
  nest_by(nsub, .drop = FALSE) %>%
  mutate(
    mean_read_n = mean(data$reads),
    sd_read_n = sd(data$reads),
    min_read_n = min(data$reads),
    max_read_n = max(data$reads)
  ) %>%
  select(nsub, mean_read_n, sd_read_n, min_read_n, max_read_n)
write_csv(summary_err_seqread_tib, paste0(path_out, "/09-summary_err_seqread_tib.csv"))

# Save R ojbects
saveRDS(sequencing, paste0(path_obj, "/sequencing.obj"))

# Tidy workspace
rm(products, read_size)

# Save workspace and session info
save.image(paste0(path_out, "/sequencing.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_out, "/sequencing.info"))
