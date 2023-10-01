# If stcpR6 is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcpR6")
library(stcpR6)
# In this analysis, we use tidyverse to handle NBA data
library(tidyverse)
# Load NBA regular season log from 2011 to 2022
dat <- read.csv("NBA_2010_2020.csv")
summary(dat$ptsTeam[dat$ptsTeam > 0])
summary(dat$plusminusTeam)
hist(dat$plusminusTeam)

# Get Cleveland Cavaliers Stats
CLE_dat <- dat %>% dplyr::filter(slugTeam == "CLE") %>%
  select(
    yearSeason,
    slugSeason,
    typeSeason,
    dateGame,
    nameTeam,
    slugTeam,
    isWin,
    ptsTeam,
    plusminusTeam
  ) %>%
  mutate(monthGame = format(as.Date(dateGame), "%Y-%m")) %>%
  mutate(relative_pm = plusminusTeam / (ptsTeam - plusminusTeam / 2))

CLE_dat_2010 <- CLE_dat %>% filter(yearSeason == 2010)
CLE_dat <-
  CLE_dat %>% filter(yearSeason > 2010 & yearSeason <= 2018)

hist(CLE_dat$plusminusTeam)


year_summ <- CLE_dat %>% group_by(yearSeason) %>%
  summarise(
    win_rate_year = mean(isWin),
    pts_year = mean(ptsTeam),
    pm_year = mean(plusminusTeam),
  )
month_summ <- CLE_dat %>% group_by(monthGame) %>%
  summarise(
    win_rate_month = mean(isWin),
    pts_month = mean(ptsTeam),
    pm_month = mean(plusminusTeam),
  )

CLE_dat <-
  CLE_dat %>% left_join(year_summ, by = "yearSeason") %>% left_join(month_summ, by = "monthGame")

regular_season_start_end_date <- CLE_dat %>%
  group_by(yearSeason) %>%
  summarise(start_date = min(dateGame),
            end_date = max(dateGame))

year_summ <-
  year_summ %>% left_join(regular_season_start_end_date, by = "yearSeason")

# NBA data ----
# 1. Win rate ----
# Pre-change : win_rate <= 0.49
# post_change: win_rate >= 0.51

plot(
  as.Date(CLE_dat$dateGame),
  CLE_dat$win_rate_month,
  pch = 20,
  xlab = "Date",
  ylab = "win rate",
  main = "Monthly win rates with seasonal averages"
)
for (i in 1:nrow(year_summ)) {
  year_date_range <-
    year_summ[i, c("start_date", "end_date")] %>% as.character() %>% as.Date()
  win_rate_year <- year_summ[i, "win_rate_year"] %>% as.numeric
  lines(
    x = year_date_range,
    y = rep(win_rate_year, 2),
    col = 2,
    lwd = 2
  )
}

# Build model
alpha <- 1e-3 # Inverse of target ARL
p_pre <- 0.49
delta_lower <- 0.02

sr_ber_model <- Stcp$new(
  method = "SR",
  family = "Ber",
  alternative = "greater",
  threshold = log(1/alpha),
  m_pre = p_pre,
  delta_lower = delta_lower,
  delta_upper = 1 - (p_pre + delta_lower)
)

cu_ber_model <- Stcp$new(
  method = "CU",
  family = "Ber",
  alternative = "greater",
  threshold = log(1/alpha),
  m_pre = p_pre,
  delta_lower = delta_lower,
  delta_upper = 1 - (p_pre + delta_lower)
)

# Compute mixtures of SR and CUSUM e-detectors.
mix_SR_ber <- sr_ber_model$updateAndReturnHistories(CLE_dat$isWin)
mix_CS_ber <- cu_ber_model$updateAndReturnHistories(CLE_dat$isWin)

# Stopping times of mixtures of e-SR and e-CUSUM procedures
mix_SR_stop_ber <- CLE_dat$dateGame[sr_ber_model$getStoppedTime()] %>% as.Date()
mix_CS_stop_ber <- CLE_dat$dateGame[cu_ber_model$getStoppedTime()] %>% as.Date()

# Plot log e-detector
# Size 600 * 450 for the paper

plot(as.Date(CLE_dat$dateGame),
     mix_SR_ber,
     xlab = "n (Date scale)",
     ylab = expression('log(M'['n']*')'),
     main = paste0("CP detected at ", mix_SR_stop_ber, " (SR) / ", mix_CS_stop_ber, " (CUSUM)"),
     type = "l",
     col = 2)
lines(as.Date(CLE_dat$dateGame),
      mix_CS_ber,
      col = 3)

# Draw customized detection line for the paper
abline(h = sr_ber_model$getThreshold())
abline(v = mix_SR_stop_ber, col = 2, lty = 2)
abline(v = mix_CS_stop_ber, col = 3, lty = 2)



# Draw detected line onto the original winning rate plot.
plot(
  as.Date(CLE_dat$dateGame),
  CLE_dat$win_rate_month,
  pch = 20,
  xlab = "Date",
  ylab = "win rate",
  main = "Monthly win rates with seasonal averages"
)
for (i in 1:nrow(year_summ)) {
  year_date_range <-
    year_summ[i, c("start_date", "end_date")] %>% as.character() %>% as.Date()
  win_rate_year <- year_summ[i, "win_rate_year"] %>% as.numeric
  lines(
    x = year_date_range,
    y = rep(win_rate_year, 2),
    col = 2,
    lwd = 2
  )
}
abline(v = mix_SR_stop_ber,
       col = 2,
       lty = 2,
       lwd = 2)




# 2. +/-  ----
# Pre-change : +/- =< -1
# Post-change: +/- > 1
# Assume +/- for each game is always between -80 and 80

plot(
  as.Date(CLE_dat$dateGame),
  CLE_dat$plusminusTeam,
  pch = 20,
  xlab = "Game Date",
  ylab = "+/-",
  main = "Plus-Minus of the Cavaliers"
)
for (i in 1:nrow(year_summ)) {
  year_date_range <-
    year_summ[i, c("start_date", "end_date")] %>% as.character() %>% as.Date()
  pm_year <- year_summ[i, "pm_year"] %>% as.numeric
  lines(
    x = year_date_range,
    y = rep(pm_year, 2),
    col = 2,
    lwd = 3
  )
}

# Build model
alpha <- 1e-3 # Inverse of target ARL
bound_lower <- -80
bound_upper <- 80
m_pre <- -1
delta_lower <- 2

apply_scale <- function(x) {
  return((x-bound_lower) / (bound_upper - bound_lower))
}

sr_bounded_model <- Stcp$new(
  method = "SR",
  family = "Bounded",
  alternative = "greater",
  threshold = log(1/alpha),
  m_pre = apply_scale(m_pre),
  delta_lower = delta_lower / (bound_upper - bound_lower)
)

cu_bounded_model <- Stcp$new(
  method = "CU",
  family = "Bounded",
  alternative = "greater",
  threshold = log(1/alpha),
  m_pre = apply_scale(m_pre),
  delta_lower = delta_lower / (bound_upper - bound_lower)
)


# Compute mixtures of SR and CUSUM e-detectors.
mix_SR_bounded <- sr_bounded_model$updateAndReturnHistories(apply_scale(CLE_dat$plusminusTeam))
mix_CS_bounded <- cu_bounded_model$updateAndReturnHistories(apply_scale(CLE_dat$plusminusTeam))

# Stopping time of the mixture of SR procedure.
mix_SR_stop_bounded <- CLE_dat$dateGame[sr_bounded_model$getStoppedTime()] %>% as.Date()
mix_CS_stop_bounded <- CLE_dat$dateGame[cu_bounded_model$getStoppedTime()] %>% as.Date()


# Plot log e-detector
# Size 600 * 450 for the paper

plot(as.Date(CLE_dat$dateGame),
     mix_SR_bounded,
     xlab = "n (Date scale)",
     ylab = expression('log(M'['n']*')'),
     main = paste0("CP detected at ", mix_SR_stop_bounded, " (SR) / ", mix_CS_stop_bounded, " (CUSUM)"),
     type = "l",
     col = 2)
lines(as.Date(CLE_dat$dateGame),
      mix_CS_bounded,
      col = 3)

# Draw customized detection line for the paper
abline(h = sr_bounded_model$getThreshold())
abline(v = mix_SR_stop_bounded, col = 2, lty = 2)
abline(v = mix_CS_stop_bounded, col = 3, lty = 2)

# M_n scale
plot(seq_along(mix_SR_bounded),
     exp(mix_SR_bounded),
     xlab = "n",
     ylab = expression('M'['n']),
     main = paste0("CP detected at ", sr_bounded_model$getStoppedTime(), "-th game (SR) / ", cu_bounded_model$getStoppedTime(), "-th game (CUSUM)"),
     type = "l",
     col = 2,
     ylim = c(0, 4 / alpha))
lines(seq_along(mix_CS_bounded),
      exp(mix_CS_bounded),
      col = 3)

# Draw customized detection line for the paper
abline(h = exp(sr_bounded_model$getThreshold()))
abline(v = sr_bounded_model$getStoppedTime(), col = 2, lty = 2)
abline(v = cu_bounded_model$getStoppedTime(), col = 3, lty = 2)



# Draw detected line onto the original +/- plot.
plot(
  as.Date(CLE_dat$dateGame),
  CLE_dat$plusminusTeam,
  pch = 20,
  xlab = "Date",
  ylab = "+/-",
  main = "+/- with detected CP"
)
for (i in 1:nrow(year_summ)) {
  year_date_range <-
    year_summ[i, c("start_date", "end_date")] %>% as.character() %>% as.Date()
  pm_year <- year_summ[i, "pm_year"] %>% as.numeric
  lines(
    x = year_date_range,
    y = rep(pm_year, 2),
    col = 2,
    lwd = 2
  )
}

abline(v = mix_SR_stop_bounded,
       col = 2,
       lty = 2,
       lwd = 2)


# Uniform mixture example
# Normalize observations
alpha <- 1e-3 # Inverse of target ARL
bound_lower <- -80
bound_upper <- 80
m_pre <- -1
pm_normalized <- apply_scale(CLE_dat$plusminusTeam)

# Build model
# Compute parameters
by_gap <- 0.02
lambda <- seq(by_gap, 1-by_gap, by = by_gap)
omega <- rep(1/length(lambda), length(lambda))

# Compute mixtures of SR and CUSUM e-detectors.
sr_bounded_uniform<- Stcp$new(
  method = "SR",
  family = "Bounded",
  alternative = "greater",
  threshold = log(1/alpha),
  m_pre = apply_scale(m_pre),
  weights = omega,
  lambdas = lambda
)

cu_bounded_uniform <- Stcp$new(
  method = "CU",
  family = "Bounded",
  alternative = "greater",
  threshold = log(1/alpha),
  m_pre = apply_scale(m_pre),
  weights = omega,
  lambdas = lambda
)

# Compute mixtures of SR and CUSUM e-detectors.
uniform_SR <- sr_bounded_uniform$updateAndReturnHistories(apply_scale(CLE_dat$plusminusTeam))
uniform_CS <- cu_bounded_uniform$updateAndReturnHistories(apply_scale(CLE_dat$plusminusTeam))

# Stopping time of the mixture of SR procedure.
uniform_SR_stop <- CLE_dat$dateGame[sr_bounded_model$getStoppedTime()] %>% as.Date()
uniform_CS_stop <- CLE_dat$dateGame[cu_bounded_model$getStoppedTime()] %>% as.Date()


plot(as.Date(CLE_dat$dateGame),
     uniform_SR,
     xlab = "n (Date scale)",
     ylab = expression('log(M'['n']*')'),
     main = paste0("CP detected at ", uniform_SR_stop, " (SR) / ", uniform_CS_stop, " (CUSUM)"),
     type = "l",
     col = 2)
lines(as.Date(CLE_dat$dateGame),
      uniform_CS,
      col = 3)
lines(as.Date(CLE_dat$dateGame),
      mix_SR_bounded,
      col = 4, lty = 2)
lines(as.Date(CLE_dat$dateGame),
      mix_CS_bounded,
      col = 5, lty = 2)

# Draw customized detection line for the paper
abline(h = sr_bounded_uniform$getThreshold())
abline(v = uniform_SR_stop, col = 2, lty = 2)
abline(v = uniform_CS_stop, col = 3, lty = 2)
abline(v = mix_SR_stop_bounded, col = 4, lty = 2)
abline(v = mix_CS_stop_bounded, col = 5, lty = 2)
