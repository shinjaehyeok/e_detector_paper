# NBA data analysis
# devtools::install_github("abresler/nbastatR")
# Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
setwd("D:/Google_drive/R/NBA_stat")
library(tidyverse)
library(nbastatR)
library(future)
plan(multiprocess)

selectedSeasons <- c(2010:2020)
T_gamelog_reg <- suppressWarnings(game_logs(seasons = selectedSeasons, league = "NBA", result_types = "team", season_types = "Regular Season"))
View(head(T_gamelog_reg))

write_csv(T_gamelog_reg, file = "NBA_2010_2020.csv")

