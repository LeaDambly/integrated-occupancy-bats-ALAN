# code to compare occupancy estimates by region
library(sparta) # Trend Analysis for Unstructured Data
library(tidyverse) # Easily Install and Load the 'Tidyverse'

# sample 999 posterior estimates of regional occupancy (averaged across grid cells),
# for all three regions of ALAN per species:year combination
path <- "D:/Lea/PhD/RProjects/ALAN/03_results/1km/"
files <- list.files(path, pattern = "*.rdata", full.names = TRUE)

modellist <- list()
for (i in 1:length(files)) {
  load(files[i])
  modellist[[i]] <- out
}

# name list item by species
for (i in 1:length(modellist)) names(modellist)[i] <- modellist[[i]]$SPP_NAME
rm(out)

getSample <- function(out, nsample = 999) {
  # This selects the posterior for the occupancy estimates (psi.fs) only
  raw_occ <- data.frame(out$BUGSoutput$sims.list)
  # Get the occupancy estimates from the summary table
  raw_occ <- raw_occ[, grep("psi.fs.r", colnames(raw_occ))]
  post <- raw_occ[sample(1:nrow(raw_occ), nsample), ]
  # Add metadata to our data
  post$spp <- out$SPP_NAME
  post$iter <- 1:nsample
  return(post)
}
post2 <- lapply(modellist, getSample)

# Function for annual growth rate
annual_growth_rate <- function(first, last, nyrs) {
  (((last / first)^(1 / nyrs)) - 1) * 100
}

# summarise by region
SumRegion <- function(post) {
  # low----
  post_low <- post %>%
    select(contains("low"))

  colnames(post_low) <- 1998:2018

  first <- post_low[, 1]
  last <- post_low[, ncol(post_low)]
  rates <- annual_growth_rate(first, last, ncol(post_low))

  ann_growth_rate <- round(mean(rates), 4)
  sd <- sd(rates)

  # Estimate the credible intervals from the posterior samples
  CI_lower <- quantile(rates, probs = 0.025, na.rm = TRUE)
  CI_upper <- quantile(rates, probs = 0.975, na.rm = TRUE)
  # Present the results
  results_low <- as.data.frame(cbind(ann_growth_rate, CI_lower, CI_upper, sd)) %>%
    mutate(region = "low") %>%
    mutate(species = post$spp[1])

  rates_low <- as.data.frame(rates) %>%
    mutate(region = "low") %>%
    mutate(species = post$spp[1])

  occupancy_low <- data.frame(occ = unlist(post_low), row.names = NULL) %>%
    mutate(region = "low") %>%
    mutate(species = post$spp[1])

  # medium----
  post_med <- post %>%
    select(contains("medium"))

  colnames(post_med) <- 1998:2018

  first <- post_med[, 1]
  last <- post_med[, ncol(post_med)]
  rates <- annual_growth_rate(first, last, ncol(post_med))

  ann_growth_rate <- round(mean(rates), 4)
  sd <- sd(rates)

  # Estimate the credible intervals from the posterior samples
  CI_lower <- quantile(rates, probs = 0.025, na.rm = TRUE)
  CI_upper <- quantile(rates, probs = 0.975, na.rm = TRUE)
  # Present the results
  results_med <- as.data.frame(cbind(ann_growth_rate, CI_lower, CI_upper, sd)) %>%
    mutate(region = "medium") %>%
    mutate(species = post$spp[1])

  rates_med <- as.data.frame(rates) %>%
    mutate(region = "medium") %>%
    mutate(species = post$spp[1])

  occupancy_med <- data.frame(occ = unlist(post_med), row.names = NULL) %>%
    mutate(region = "medium") %>%
    mutate(species = post$spp[1])

  # high----
  post_high <- post %>%
    select(contains("high"))

  colnames(post_high) <- 1998:2018

  first <- post_high[, 1]
  last <- post_high[, ncol(post_high)]
  rates <- annual_growth_rate(first, last, ncol(post_high))

  ann_growth_rate <- round(mean(rates), 4)
  sd <- sd(rates)

  # Estimate the credible intervals from the posterior samples
  CI_lower <- quantile(rates, probs = 0.025, na.rm = TRUE)
  CI_upper <- quantile(rates, probs = 0.975, na.rm = TRUE)
  # Present the results
  results_high <- as.data.frame(cbind(ann_growth_rate, CI_lower, CI_upper, sd)) %>%
    mutate(region = "high") %>%
    mutate(species = post$spp[1])

  rates_high <- as.data.frame(rates) %>%
    mutate(region = "high") %>%
    mutate(species = post$spp[1])

  occupancy_high <- data.frame(occ = unlist(post_high), row.names = NULL) %>%
    mutate(region = "high") %>%
    mutate(species = post$spp[1])

  # all together now----
  results <- bind_rows(results_low, results_med, results_high)
  rownames(results) <- NULL
  rates <- bind_rows(rates_low, rates_med, rates_high)
  occupancy <- bind_rows(occupancy_low, occupancy_med, occupancy_high)
  ll <- list(results, rates, occupancy)
  return(ll)
}

ll <- lapply(post2, SumRegion)

# Does occupancy differ in ALAN regions within species?----
# serotinus
e.ser <- ll[[1]][[3]]

# is Medium larger than Low?
dif1 <- e.ser %>%
  filter(region != "high") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = medium - low)

# get the proportion of values in diff that are positive
# value from 0 to 1
# which can be interpreted as a confidence that A is greater than B
dif1 <- dif1$diff
round(mean(dif1 > 0), 2)

# is Medium larger than High?
dif2 <- e.ser %>%
  filter(region != "low") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = medium - high)

dif2 <- dif2$diff
round(mean(dif2 > 0), 2)

# is Low larger than High?
dif3 <- e.ser %>%
  filter(region != "medium") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = low - high)

dif3 <- dif3$diff
round(mean(dif3 > 0), 2)

# noctule
n.noc <- ll[[2]][[3]]

# is High larger than Medium?
dif4 <- n.noc %>%
  filter(region != "low") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = high - medium)

dif4 <- dif4$diff
round(mean(dif4 > 0), 2)

# is High larger than Low?
dif5 <- n.noc %>%
  filter(region != "medium") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = high - low)

dif5 <- dif5$diff
round(mean(dif5 > 0), 2)

# is Medium larger than Low?
dif6 <- n.noc %>%
  filter(region != "high") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = medium - low)

dif6 <- dif6$diff
round(mean(dif6 > 0), 2)

# pip1
p.pip <- ll[[3]][[3]]

# is High larger than Medium?
dif7 <- p.pip %>%
  filter(region != "low") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = high - medium)

dif7 <- dif7$diff
round(mean(dif7 > 0), 2)

# is High larger than Low?
dif8 <- p.pip %>%
  filter(region != "medium") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = high - low)

dif8 <- dif8$diff
round(mean(dif8 > 0), 2)

# is Medium larger than Low?
dif9 <- p.pip %>%
  filter(region != "high") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = medium - low)

dif9 <- dif9$diff
round(mean(dif9 > 0), 2)

# pip2
p.pyg <- ll[[4]][[3]]

# is High larger than Medium?
dif10 <- p.pyg %>%
  filter(region != "low") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = high - medium)

dif10 <- dif10$diff
round(mean(dif10 > 0), 2)

# is High larger than Low?
dif11 <- p.pyg %>%
  filter(region != "medium") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = high - low)

dif11 <- dif11$diff
round(mean(dif11 > 0), 2)

# is Medium larger than Low?
dif12 <- p.pyg %>%
  filter(region != "high") %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = "region", values_from = "occ") %>%
  select(-id) %>%
  mutate(diff = medium - low)

dif12 <- dif12$diff
round(mean(dif12 > 0), 2)
