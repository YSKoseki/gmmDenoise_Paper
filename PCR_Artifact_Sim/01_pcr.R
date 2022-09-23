# 01_pcr.R
# A simulation of PCR amplification of true and error sequences

# Increase vector memory limit: Save below as the first line of '.Renviron'
#   'R_MAX_VSIZE = 200Gb'
#   Note that the value required may vary by machine
#usethis::edit_r_environ("project")

# Packages required
library(dirmult); packageVersion("dirmult"); # 0.1.3.5
library(doParallel); packageVersion("doParallel") # 1.0.14
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Parallel computation setup
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

# Load helper functions
source("./99_helper.R")

# Parameters 
repl_n <- 100
ini_seq_n <- 100
ini_seq_len <- 160
ini_seq_GC <- 56
ini_read_n <- 10000
dirich_alp <- 1
cyc <- 35
beta_alp1 <- beta_beta1 <- 5 # See Kelly et al., 2019
meanlog_e1 <- 0 # A small amount of stochasticity (Kelly et al., 2019)
sdlog_e1 <- .05 # The same as above
mean_e2 <- 1.2e-05 # Following the estimates for KOD in Potapov & Ong (2017)
sd_e2 <- .2e-05 # The same as above
beta_params(mean = mean_e2, sd = sd_e2) # 3.599956e+01 2.999927e+06
beta_alp2 <- 36
beta_beta2 <- 3e+06

# Output directories and files
path_out <- "./01-pcr.R"
dir.create(path_out, recursive = TRUE)
path_repl <- paste0(path_out, "/repl_",
                     formatC(1:repl_n, width = nchar(repl_n), flag = "0"))
path_repl %>% lapply(dir.create, recursive = TRUE)
path_prod <- 0:cyc %>% formatC(width = nchar(cyc), flag = "0") %>%
  paste0("/products_", ., ".csv") %>%
  expand_grid(folder = path_repl, file = .) %>%
  mutate(path = paste0(folder, file)) %>% select(path) %>% pull() %>%
  matrix(cyc + 1, repl_n)
path_read <- paste0(path_repl, "/read_size.csv")

# Writing initial data set into csv
foreach (j = 1:repl_n, .packages = c("dirmult", "tidyverse")) %dopar% {
  set.seed(j)
  products_00 <- read_size <- NULL
  # Template data
  products_00 <- 0 %>%
    expand_grid(
      repl = j,
      cyc = .,
      seq = 1:ini_seq_n %>%
        formatC(width = nchar(ini_seq_n), flag = "0") %>%
        paste0("Seq", .)
    ) %>%
    mutate(
      conc = dirmult::rdirichlet(
        n = 1, alpha = rep(dirich_alp, ini_seq_n)
      ) %>% t() %>% c(),
      # Amplification efficiency
      eff = rbeta(
        n = ini_seq_n, shape1 = beta_alp1, shape2 = beta_beta1
      ),
      dna = random_seq(
        n = ini_seq_n, seq_len = ini_seq_len, GC_ratio = ini_seq_GC
      ) %>% as.character()
    )
  # Read size data
  read_size <- tibble(
    repl = j,
    cyc = 0,
    N = ini_read_n,
    sum_conc = sum(products_00$conc),
    sum_true_conc = 1,
    sum_err_conc = 0
  )
  write_csv(products_00, path_prod[1, j])
  write_csv(read_size, path_read[j])
}
gc(verbose = TRUE, reset = TRUE)

# PCR Simulation
foreach (j = 21, .packages = "tidyverse") %do% {
  set.seed(j)
  for (i in 1:cyc) {
    products <- read_size <- sum_conc <- N <- amplicons <- ampl_n <-
      error_pos <- errorseq_char <- errorseq_nam <- errorseq_conc <-
      errorseq_eff <- error_pos_ul <- flag <- errorseq_n <- errorseqs <-
      correctseqs <- sum_true_conc <- sum_err_conc <- NULL
    # Read data set
    products <- read_csv(path_prod[i, j], show_col_types = FALSE)
    read_size <- read_csv(path_read[j], show_col_types = FALSE)
    # Amplify sequences at a beta-distributed efficiency with a log-normal error
    sum_conc <- read_size %>%
      filter(cyc == (i - 1)) %>% select(sum_conc) %>% pull()
    N <- read_size %>%
      filter(cyc == (i - 1)) %>% select(N) %>% pull()
    amplicons <- products %>%
      mutate(e1 = nrow(products) %>% rlnorm(meanlog_e1, sdlog_e1),
             conc_ampl = (conc / sum_conc) * eff * e1,
             N_ampl = (N * conc_ampl) %>% floor())
    # Generate amplification error (base substitution) at a beta-distributed rate
    ampl_n <- nrow(amplicons)
    amplicons <- amplicons %>%
      mutate(e2 = rbeta(ampl_n, shape1 = beta_alp2, shape2 = beta_beta2),
             conc_err = conc_ampl * ini_seq_len * e2,
             N_err = (N * conc_err) %>% floor())
    # Choose base positions of substitution, assuming that substitution occurs
    #   only once in a sequence given the orders of error rates
    error_pos <- amplicons$N_err %>%
      lapply(function(x) {
        if (x > 0) {
          y <- rdunif(x, b = ini_seq_len, a = 1)
        } else {
          y <- NULL
        }
        return(y)
      })
    # Substitute the assigned bases, assuming that only transition occurs
    #   (cf. McInerney et al., 2014; Potapov & Ong, 2017)
    errorseq_char <- amplicons$dna %>% as.list() %>% 
      map2(error_pos, function(x, y) {
        if (!is.null(y)) {
          z <- transition(x, y)
        } else {
          z <- NULL
        }
        return(z)
      })
    # Name the newly generated error sequences
    errorseq_nam <- amplicons$seq %>% as.list() %>%
      map2(error_pos, function(x, y) {
        if (!is.null(y)) {
          n <- length(y)
          z <- paste0(x, "e") %>% rep(n)
        } else {
          z <- NULL
        }
        return(z)
      })
    # Split error concentration evenly among homologous errors
    errorseq_conc <- amplicons$conc_err %>% as.list() %>%
      map2(error_pos, function(x, y) {
        if (!is.null(y)) {
          n <- length(y)
          z <- x / n
          v <- rep(z, n)
        } else {
          v <- NULL
        }
        return(v)
      })
    # Error sequences are assumed to have the same values of amplification
    #   efficiency to those of the templates
    errorseq_eff <- amplicons$eff %>% as.list() %>%
      map2(error_pos, function(x, y) {
        if (!is.null(y)) {
          n <- length(y)
          z <- rep(x, n)
        } else {
          z <- NULL
        }
        return(z)
      })
    # Tabulate error sequence data
    error_pos_ul <- unlist(error_pos)
    flag <- !is.null(error_pos_ul)
    errorseq_n <- length(error_pos_ul)
    errorseqs <- if (flag) {
      expand_grid(repl = first(amplicons$repl),
                  cyc = first(amplicons$cyc),
                  seq = unlist(errorseq_nam)) %>%
        mutate(N = sum(amplicons$N_ampl),
               conc = unlist(errorseq_conc),
               eff = unlist(errorseq_eff),
               dna = unlist(errorseq_char))
    } else {
      NULL
    }
    # Tabulate correct sequence data
    correctseqs <- amplicons %>%
      mutate(conc = if_else(N_err > 0,
                            conc + conc_ampl - conc_err,
                            conc + conc_ampl)) %>%
      select(repl, cyc, seq, conc, eff, dna)
    # Sort amplicon sequences (correct and error) and write the
    #   sorted sequence data into csv
    products <- NULL
    products <- bind_rows(correctseqs, errorseqs) %>%
      group_by(dna) %>%
      summarize(repl = j, cyc = i, seq = first(seq), conc = sum(conc),
                eff = first(eff), dna= first(dna)) %>%
      select(repl, cyc, seq, conc, eff, dna) %>%
      arrange(seq, desc(conc))
    write_csv(products, path_prod[i + 1, j])
    # Update the total read size and concentration data in csv
    sum_true_conc <- products %>%
      filter(!str_detect(.$seq, "e$")) %>% select(conc) %>% pull() %>% sum()
    sum_err_conc <- products %>%
      filter(str_detect(.$seq, "e$")) %>% select(conc) %>% pull() %>% sum()
    sum(amplicons$conc_ampl)
    read_size <- tibble(repl = j,
                        cyc = i,
                        N = N + sum(amplicons$N_ampl),
                        sum_conc = sum_conc + sum(amplicons$conc_ampl),
                        sum_true_conc = sum_true_conc,
                        sum_err_conc = sum_err_conc) %>%
      bind_rows(read_size, .)
    write_csv(read_size, path_read[j])
    # Garbage collection
    rm(products, read_size, sum_conc, N, amplicons, ampl_n, error_pos,
       errorseq_char, errorseq_nam, errorseq_conc, errorseq_eff, error_pos_ul,
       flag, errorseq_n, errorseqs, correctseqs, sum_true_conc, sum_err_conc)
    gc(verbose = TRUE, reset = TRUE)
  }
}
rm(i, j)
# End of parallel computation
stopCluster(cl)

# Save workspace and session info
save.image(paste0(path_out, "/pcr.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_out, "/pcr.info"))


