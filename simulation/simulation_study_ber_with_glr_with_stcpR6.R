# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcpR6")

library(stcpR6)
# Bernoulli case ----
p_pre <- 0.5
p_post <- 0.6
min_delta <- 0.01
max_sample <- 1000L
v <- 0 # Worst-case for SR type detector
ARL_target <- max_sample * 0.5
alpha <- 1 / ARL_target
num_repeat <- 5000L

# CUSUM functions
run_ber_cusum_with_reset_thres <- function(x_vec, thres) {
  oracle_cu <- Stcp$new(
    method = "CU",
    family = "Ber",
    alternative = "greater",
    threshold = thres,
    m_pre = p_pre,
    delta_lower = p_post - p_pre,
    delta_upper = p_post - p_pre
  )
  oracle_cu$updateLogValuesUntilStop(x_vec)
  return(oracle_cu$getStoppedTime())
}

# GLR CUSUM
run_ber_glr_cusum_with_reset_thres <- function(x_vec, thres) {
  glr_cu <- Stcp$new(
    method = "GLRCU",
    family = "Ber",
    alternative = "greater",
    threshold = thres,
    m_pre = p_pre
  )
  glr_cu$updateLogValuesUntilStop(x_vec)
  return(glr_cu$getStoppedTime())
}

# E-detectors
run_ber_eCU_with_reset_thres <- function(x_vec, thres) {
  stcp_cu <- Stcp$new(
    method = "CU",
    family = "Ber",
    alternative = "greater",
    threshold = thres,
    m_pre = p_pre,
    delta_lower = min_delta,
    delta_upper = 1 - min_delta - p_pre
  )
  stcp_cu$updateLogValuesUntilStop(x_vec)
  return(stcp_cu$getStoppedTime())
}


# Simulation functions
AD_ber_fn <-
  function(thres, v, run_fn) {
    single_run <- function() {
      x_vec <- c(rbinom(v, 1, p_pre), rbinom(max_sample - v, 1, p_post))
      stopped_time <- run_fn(x_vec, thres)
      run_length <-
        ifelse(stopped_time == 0, max_sample, stopped_time)
      return(run_length)
    }
    simulation_with_fixed_thres <-
      replicate(num_repeat, single_run())
    return(list(
      v = v,
      m = mean(simulation_with_fixed_thres[simulation_with_fixed_thres > v]) - v,
      sd = sd(simulation_with_fixed_thres[simulation_with_fixed_thres > v]),
      simulation_raw = simulation_with_fixed_thres
    ))
  }

WAD_ber_fn <- function(thres,
                       run_fn,
                       v_vec = seq(0, max_sample / 2, length.out = 6L)) {
  f <- function(v) {
    AD_ber_fn(thres, v, run_fn)
  }
  return(sapply(v_vec, f))
}

ARL_ber_fn <- function(thres,
                       run_fn) {
  single_run <- function() {
    x_vec <- rbinom(max_sample, 1, p_pre)
    stopped_time <- run_fn(x_vec, thres)
    run_length <-
      ifelse(stopped_time == 0, max_sample, stopped_time)
    return(run_length)
  }
  simulation_with_fixed_thres <- replicate(num_repeat, single_run())
  return(list(
    m = mean(simulation_with_fixed_thres),
    sd = sd(simulation_with_fixed_thres)
  ))
}

find_exact_thres <- function(run_fn) {
  f <- function(thres) {
    arl <-
      ARL_ber_fn(thres, run_fn)
    return(round(arl$m - ARL_target - 1))
  }
  exact_thres <- stats::uniroot(f,
                                c(0.1, log(ARL_target) * 10), tol = 1e-3, maxiter = 100)
  return(exact_thres$root)
}

# 1. Exact CUSUM.
format(Sys.time(), "%c")
set.seed(100)
exact_cusum_thres <- find_exact_thres(run_fn = run_ber_cusum_with_reset_thres)

oracle_cu <- Stcp$new(
  method = "CU",
  family = "Ber",
  alternative = "greater",
  threshold = exact_cusum_thres,
  m_pre = p_pre,
  delta_lower = p_post - p_pre,
  delta_upper = p_post - p_pre
)

run_ber_cusum <- function(x_vec, ...) {
  oracle_cu$reset()
  oracle_cu$updateLogValuesUntilStop(x_vec)
  return(oracle_cu$getStoppedTime())
}

ARL_exact_cusum <- ARL_ber_fn(
  thres = exact_cusum_thres,
  run_fn = run_ber_cusum
)

WAD_exact_cusum <- WAD_ber_fn(
  thres = exact_cusum_thres,
  run_fn = run_ber_cusum
)

# 2. Exact GLR CUSUM.
format(Sys.time(), "%c")
set.seed(100)
exact_glr_cusum_thres <- find_exact_thres(run_fn = run_ber_glr_cusum_with_reset_thres)

glr_cu <- Stcp$new(
  method = "GLRCU",
  family = "Ber",
  alternative = "greater",
  threshold = exact_glr_cusum_thres,
  m_pre = p_pre
)

run_ber_glr_cusum <- function(x_vec, ...) {
  glr_cu$reset()
  glr_cu$updateLogValuesUntilStop(x_vec)
  return(glr_cu$getStoppedTime())
}

ARL_exact_glr_cusum <- ARL_ber_fn(
  thres = exact_glr_cusum_thres,
  run_fn = run_ber_glr_cusum
)
WAD_exact_glr_cusum <- WAD_ber_fn(
  thres = exact_glr_cusum_thres,
  run_fn = run_ber_glr_cusum
)

# 3. e-SR with unknown p_post
format(Sys.time(), "%c")
set.seed(100)

stcp_sr <- Stcp$new(
  method = "SR",
  family = "Ber",
  alternative = "greater",
  threshold = log(ARL_target),
  m_pre = p_pre,
  delta_lower = min_delta,
  delta_upper = 1 - min_delta - p_pre
)

run_ber_eSR <- function(x_vec,...) {
  stcp_sr$reset()
  stcp_sr$updateLogValuesUntilStop(x_vec)
  return(stcp_sr$getStoppedTime())
}

ARL_eSR_mix <- ARL_ber_fn(
  thres = log(ARL_target),
  run_fn = run_ber_eSR
)
WAD_eSR_mix <- WAD_ber_fn(
  thres = log(ARL_target),
  run_fn = run_ber_eSR
)

# 4. e-CUSUM with unknown p_post and optimal threshold
format(Sys.time(), "%c")
set.seed(100)
exact_eCS_thres <- find_exact_thres(run_fn = run_ber_eCU_with_reset_thres)

stcp_cu <- Stcp$new(
  method = "CU",
  family = "Ber",
  alternative = "greater",
  threshold = exact_eCS_thres,
  m_pre = p_pre,
  delta_lower = min_delta,
  delta_upper = 1 - min_delta - p_pre
)

run_ber_eCU <- function(x_vec, ...) {
  stcp_cu$reset()
  stcp_cu$updateLogValuesUntilStop(x_vec)
  return(stcp_cu$getStoppedTime())
}

ARL_exact_eCS_mix <- ARL_ber_fn(
  thres = exact_eCS_thres,
  run_fn = run_ber_eCU
)
WAD_exact_eCS_mix <- WAD_ber_fn(
  thres = exact_eCS_thres,
  run_fn = run_ber_eCU
)


# Single run visualization

single_run_plot <- function(x_vec, v = NA, plot_raw = FALSE) {
  
  oracle_cu$reset()
  glr_cu$reset()
  stcp_sr$reset()
  stcp_cu$reset()
  
  exact_cusum <- oracle_cu$updateAndReturnHistories(x_vec)
  exact_glr_cusum <- glr_cu$updateAndReturnHistories(x_vec)

  eCS_mix <- stcp_sr$updateAndReturnHistories(x_vec)
  eSR_mix <- stcp_cu$updateAndReturnHistories(x_vec)

  v_inner <- ifelse(is.na(v), Inf, v)

  if (plot_raw) {
    plot(
      seq_along(x_vec),
      exact_cusum,
      type = "l",
      ylim = c(0, 3 * stcp_sr$getThreshold()),
      xlab = "n",
      ylab = "CP-stat"
    )
    abline(v = v_inner, lty = 2)
    abline(h = exact_cusum_thres)
    lines(seq_along(x_vec), exact_glr_cusum, col = 2)
    abline(h = exact_glr_cusum_thres, col = 2)
    lines(seq_along(x_vec), eCS_mix, col = 3)
    abline(h = exact_eCS_thres, col = 3)
    lines(seq_along(x_vec), eSR_mix, col = 4)
    abline(h = stcp_sr$getThreshold(),
           col = 4)
  } else {
    plot(
      seq_along(x_vec),
      exact_cusum / exact_cusum_thres,
      type = "l",
      ylim = c(0, 3),
      xlab = "n",
      ylab = "CP-stat"
    )
    abline(v = v_inner, lty = 2)
    abline(h = 1)
    lines(seq_along(x_vec),
          exact_glr_cusum / exact_glr_cusum_thres,
          col = 2)
    lines(seq_along(x_vec),
          eCS_mix / exact_eCS_thres,
          col = 3)
    lines(
      seq_along(x_vec),
      eSR_mix / stcp_sr$getThreshold(),
      col = 4
    )
  }
}
# No change
set.seed(100)
x_vec <- rbinom(ARL_target, 1, p_pre)
single_run_plot(x_vec)

# v = 0
set.seed(100)
v <- 0
x_vec <- c(rbinom(v, 1, p_pre), rbinom(ARL_target - v, 1, p_post))
single_run_plot(x_vec, v)

# v = 200
set.seed(100)
v <- 200
x_vec <- c(rbinom(v, 1, p_pre), rbinom(ARL_target - v, 1, p_post))
single_run_plot(x_vec, v)
single_run_plot(x_vec, v, plot_raw = TRUE)
