library(tidyverse)

num_repeat <- 5000L

get_simul_tb <- function(simul_list, method) {
  get_each_run_tb <- function(i) {
    data.frame(
      v = simul_list["v",i][[1]],
      stopped_time = simul_list["simulation_raw",i][[1]]
    )
  }
  simul_raw_tb <- do.call("rbind", lapply(1:dim(simul_list)[2], get_each_run_tb)) |>
    dplyr::mutate(
      delay = stopped_time - v
    )
  summ <- simul_raw_tb |>
    dplyr::group_by(v) |>
    dplyr::summarise(
      AD = mean(delay[delay>0]),
      sd = sd(delay[delay>0]),
      ci_width = sd * 1.96 / sqrt(num_repeat),
      early_stop_ratio = 1-mean(delay>0),
      Method = method
    ) |>
    dplyr::ungroup()

  simul_raw_tb$v <- as.factor(simul_raw_tb$v)

  return(list(
    raw_tb = simul_raw_tb,
    summ = summ
  ))
}
oracle_cusum_tb <- get_simul_tb(WAD_exact_cusum, "1. Oracle CUSUM with exact threshold")
glr_cusum_tb <- get_simul_tb(WAD_exact_glr_cusum, "2. GLR with exact threshold")
exact_ecs_tb <- get_simul_tb(WAD_exact_eCS_mix, "3. e-CUSUM with exact threshold")
esr_tb <- get_simul_tb(WAD_eSR_mix, "4. e-SR with log(ARL) threshold")

summ_tb_all <- oracle_cusum_tb$summ |>
  rbind(glr_cusum_tb$summ) |>
  rbind(esr_tb$summ) |>
  rbind(exact_ecs_tb$summ)

# 900 * 500 for paper 
summ_tb_all |>
  ggplot(aes(x=v, y=AD, color = Method, shape = Method, linetype=Method)) +
    geom_errorbar(aes(ymin=AD-ci_width, ymax=AD+ci_width), width=.1) +
    geom_line() +
    geom_point(size = 2.5) +
    ylab("Average detection delay") +
    theme(
      axis.title = element_text(size = 12),
      legend.title = element_text(size=12),
      legend.text = element_text(size=10)
      ) 

# 900 * 500 for paper 
summ_tb_all |>
  ggplot(aes(x=v, y=early_stop_ratio, color = Method, shape = Method, linetype=Method)) +
  geom_line() +
  geom_point(size = 2.5) +
  ylab("Pre-change false alarm rate") +
  theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  ) 

summ_tb_all
