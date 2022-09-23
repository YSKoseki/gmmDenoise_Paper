# 99_helper.R
# A set of helper functions for simulating PCR processes

# Generate a set of random sequences
random_seq <- function(n = 10, seq_len = 100, GC_ratio = 56,
                       id_prefix = "Seq_") {
  # DNA pool under the GC_ratio 
  bases <- c("A", "T", "G", "C")
  AT_ratio <- 100 - GC_ratio
  acgt_ratio <- c(AT_ratio, AT_ratio, GC_ratio, GC_ratio) / 2
  pool <- rep(bases, acgt_ratio)
  # Generate random sequences
  delta <- n
  while (delta > 0) {
    seq <- rep(NA, n)
    seq <- sample(pool, seq_len * n, replace = TRUE) %>%
      matrix(c(seq_len, n)) %>%
      apply(2, function(x) {paste(x, collapse = "")})
    uniq_seq <- seq %>% unique()
    n_uniq <- uniq_seq %>% length()
    delta <- n - n_uniq
  }
  seq <- seq %>% Biostrings::DNAStringSet()
  # # Sequence IDs
  # id <- paste0(id_prefix, formatC(1:n, width = nchar(n), flag = "0"))
  # names(seq) <- id
  return(seq)
}

# Single-base transition
transition <- function(DNA, pos) {
  pre <- str_sub(DNA, start = pos, end = pos)
  aft <- case_when(pre == "A" ~ "G", pre == "G" ~ "A",
                   pre == "C" ~ "T", pre == "T" ~ "C")
  DNAaft <- DNA
  str_sub(DNAaft, start = pos, end = pos) <- aft
  return(DNAaft)
}

# Generate log-normal random numbers specifying mean and sd of the log-normal
#   distribution, not the distribution on the log scale
rlnorm2 <- function(n, mean = 1, sd = 1) {
  sdlog <- sqrt(log((sd / mean) ^ 2 + 1))
  meanlog <- log(mean) - (sdlog ^ 2) / 2
  val <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  return(val)
}

# Calculate mean of log-normal distribution whose logarithm has mean
#   equal to `meanlog` and SD equal to `sdlog`
lognorm_mean <- function(meanlog, sdlog) {
  m <- exp(meanlog + (sdlog ^ 2) / 2)
  return(m)
}

# Calculate SD of log-normal distribution whose logarithm has mean
#   equal to `meanlog` and SD equal to `sdlog`
lognorm_sd <- function(meanlog, sdlog) {
  sd <- sqrt(exp(2 * meanlog + (sdlog ^ 2)) * (exp(sdlog ^ 2) - 1))
  return(sd)
}

# Calculate log-scale mean of log-normal distribution whose mean
#   and SD are equal to `mean` and `sd`, respectively
logscale_mean <- function(mean, sd) {
  sigma <- sqrt(log((sd / mean) ^ 2 + 1))
  mu <- log(mean) - (sigma ^ 2) / 2
  return(mu)
}

# Calculate log-scale SD of log-normal distribution whose mean
#   and SD are equal to `mean` and `sd`, respectively
logscale_sd <- function(mean, sd) {
  sigma <- sqrt(log((sd / mean) ^ 2 + 1))
  return(sigma)
}

# Estimate beta distribution parameters from expected value (mean) and variance
beta_params <- function(mean, sd) {
  alpha <- ((mean ^ 2) * (1 - mean) / (sd ^ 2)) - mean
  beta <- alpha / mean - alpha
  params <- c(alpha, beta)
  names(params) <- c("alpha", "beta")
  return(params)
}
