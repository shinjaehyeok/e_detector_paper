# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")

library(stcp)
# Bernoulli case ----
p_pre <- 0.5
p_post <- 0.6
min_delta <- 0.01
max_sample <- 1000L
v <- 0 # Worst-case for SR type detector
ARL_target <- max_sample * 0.5
alpha <- 1 / ARL_target
num_repeat <- 500L

# CUSUM functions
run_ber_cusum <- function(x_vec, p_pre, p_post, thres) {
  n_max <- length(x_vec)
  n_star <- Inf
  m <- 0
  for (i in 1:n_max) {
    x <- x_vec[i]
    # Update the change statistic
    m <- max(m, 0) + ifelse(x == 1, log(p_post/p_pre), log((1-p_post)/(1-p_pre)))
    # Check whether the stopping time happens
    if (m > thres) {
      n_star <- i
      break
    }
  }
  return(n_star)
}

run_ber_cusum_with_memory <- function(x_vec, p_pre, p_post, thres) {
  n_max <- length(x_vec)
  n_star <- Inf
  m <- 0
  m_vec <- length(n_max)
  for (i in 1:n_max) {
    x <- x_vec[i]
    # Update the change statistic
    m <- max(m, 0) + ifelse(x == 1, log(p_post/p_pre), log((1-p_post)/(1-p_pre)))
    m_vec[i] <- m
    # Check whether the stopping time happens
    if (m > thres) {
      n_star <- i
    }
  }
  return(list(n_star = n_star, m_vec = m_vec))
}

# GLR CUSUM
run_ber_glr_cusum <- function(x_vec, p_pre, p_post, thres, p_min_delta = min_delta, window_size = length(x_vec)) {
  n_max <- length(x_vec)
  n_star <- Inf
  p_lower <- p_pre + p_min_delta
  m <- -Inf
  for (i in 1:n_max) {
    for (j in 1:i) {
      n_inner <- i - j + 1
      if (n_inner > window_size) {
        next
      }
      s <- sum(x_vec[j:i])
      f <- n_inner - s
      p_hat <- min(max(s/n_inner, p_lower), 0.99)
      m_inner <-  s * log(p_hat/p_pre) + f * log((1-p_hat)/(1-p_pre))
      if (m < m_inner) {
        m <- m_inner
      }
    }
    if (m > thres) {
      n_star <- i
      break
    }
    m <- -Inf
  }
  return(n_star)
}

run_ber_glr_cusum_with_memory <- function(x_vec, p_pre, p_post, thres, p_min_delta = min_delta, window_size = length(x_vec)) {
  n_max <- length(x_vec)
  n_star <- Inf
  m <- -Inf
  p_lower <- p_pre + p_min_delta
  p_hat <- p_lower
  m_vec <- length(n_max)
  p_hat_vec <- length(n_max)
  for (i in 1:n_max) {
    for (j in 1:i) {
      n_inner <- i - j + 1
      s <- sum(x_vec[j:i])
      f <- n_inner - s
      p_hat_inner <- min(max(s/n_inner, p_lower), 0.99)
      m_inner <-  s * log(p_hat_inner/p_pre) + f * log((1-p_hat_inner)/(1-p_pre))
      # print("i")
      # print(i)
      # print("j")
      # print(j)
      # print("m")
      # print(m_inner)
      # print("phat")
      # print(p_hat_inner)
      if (m < m_inner) {
        m <- m_inner
        p_hat <- p_hat_inner
      }
    }
    m_vec[i] <- m
    p_hat_vec[i] <- p_hat
    if (m > thres) {
      n_star <- i
    }
    m <- -Inf
  }
  return(list(n_star = n_star, m_vec = m_vec, p_hat_vec = p_hat_vec))
}

# E-detectors
# Compute optimal delta star
delta_star <- (p_post - p_pre)
delta_upper <- 1 - min_delta - p_pre
delta_lower <- min_delta

# Other settings
psi_fn_list <- generate_sub_B_fn(p = p_pre)
v_min <- 1
k_max <- 1e+3

# Build CP detectors
# When delta_lower = delta_upper = delta_star
stcp_star <- build_stcp_exp(
  alpha,
  p_pre,
  delta_star,
  delta_star,
  is_test = FALSE,
  psi_fn_list,
  s_fn = function(x) {
    x - p_pre
  },
  v_fn = function(x) {
    1
  },
  v_min,
  k_max,
  tol = 1e-6
)

# When delta_lower < delta_star < delta_upper
stcp_mix <- build_stcp_exp(
  alpha,
  p_pre,
  delta_lower,
  delta_upper,
  is_test = FALSE,
  psi_fn_list,
  s_fn = function(x) {
    x - p_pre
  },
  v_fn = function(x) {
    1
  },
  v_min,
  k_max,
  tol = 1e-6
)

fn_factor_run_ber_e_detector <- function(stcp_obj, is_SR_type) {
  run_ber_e_detector <- function(x_vec, p_pre, p_post, thres) {
    n_max <- length(x_vec)
    n_star <- Inf
    m <- 0
    for (i in 1:n_max) {
      x <- x_vec[i]
      # Update the change statistic
      m <- update_log_mix_e(
        x,
        stcp_obj$omega,
        stcp_obj$log_base_fn_list,
        stcp_obj$log_e_vec,
        stcp_obj$is_test,
        is_SR_type
      )
      stcp_obj$log_e_vec <- m$last_log_e_vec
      stcp_obj$n <- stcp_obj$n + 1
      # Check whether the stopping time happens
      if (m$log_mix_e_vec > thres) {
        n_star <- i
        break
      }
    }
    return(n_star)
  }
  return(run_ber_e_detector)
}

run_ber_eSR_mix <- fn_factor_run_ber_e_detector(stcp_mix, is_SR_type = TRUE)
run_ber_eCS_mix <- fn_factor_run_ber_e_detector(stcp_mix, is_SR_type = FALSE)


# Simulation functions
AD_ber_fn <- function(p_pre, p_post, thres, v, max_sample, num_repeat,
                          run_fn = run_ber_cusum) {
  single_run <- function() {
    x_vec <- c(rbinom(v, 1, p_pre), rbinom(max_sample - v, 1, p_post))
    stopped_time <- run_fn(x_vec, p_pre, p_post, thres)
    run_length <- ifelse(is.infinite(stopped_time), max_sample, stopped_time)
    return(run_length)
  }
  simulation_with_fixed_thres <- replicate(num_repeat, single_run())
  return(
    list(
      v = v,
      m = mean(simulation_with_fixed_thres[simulation_with_fixed_thres > v]) - v,
      sd = sd(simulation_with_fixed_thres[simulation_with_fixed_thres > v]),
      simulation_raw = simulation_with_fixed_thres
      )
    )
}

WAD_ber_fn <- function(p_pre, p_post,thres, max_sample, num_repeat,
                       run_fn = run_ber_cusum,
                       v_vec = seq(0, max_sample/2, length.out = 6L)) {
  f <- function(v) {
    AD_ber_fn(p_pre, p_post, thres, v, max_sample, num_repeat, run_fn)
  }
  return(sapply(v_vec, f))
}

ARL_ber_fn <- function(p_pre, p_post, thres, max_sample, num_repeat,
                          run_fn = run_ber_cusum) {
  single_run <- function() {
    x_vec <- rbinom(max_sample, 1, p_pre)
    stopped_time <- run_fn(x_vec, p_pre, p_post, thres)
    run_length <- ifelse(is.infinite(stopped_time), max_sample, stopped_time)
    return(run_length)
  }
  simulation_with_fixed_thres <- replicate(num_repeat, single_run())
  return(
    list(
      m = mean(simulation_with_fixed_thres),
      sd = sd(simulation_with_fixed_thres)
    )
  )
}

find_exact_thres <- function(p_pre, p_post, ARL_target, max_sample, num_repeat,
                             run_fn = run_ber_cusum) {
  f <- function(thres) {
    arl <- ARL_ber_fn(p_pre, p_post, thres, max_sample, num_repeat, run_fn)
    return(round(arl$m - ARL_target - 1))
  }
  exact_thres <- stats::uniroot(f,
                             c(0.1, log(ARL_target) * 10), tol = 1e-3, maxiter = 100)
  return(exact_thres$root)
}

# 1. Exact CUSUM.
format(Sys.time(), "%c")
set.seed(100)
exact_cusum_thres <- find_exact_thres(p_pre, p_post,
                                      ARL_target, max_sample, num_repeat)
ARL_exact_cusum <- ARL_ber_fn(p_pre, p_post,
                                 thres = exact_cusum_thres,
                                 max_sample = max_sample,
                                 num_repeat = num_repeat)
WAD_exact_cusum <- WAD_ber_fn(p_pre, p_post,
                                 thres = exact_cusum_thres,
                                 max_sample = max_sample,
                                 num_repeat = num_repeat)

# 2. Exact GLR CUSUM.
format(Sys.time(), "%c")
set.seed(100)
exact_glr_cusum_thres <- find_exact_thres(p_pre, p_post,
                                      ARL_target, max_sample, num_repeat,
                                      run_fn = run_ber_glr_cusum)
ARL_exact_glr_cusum <- ARL_ber_fn(p_pre, p_post,
                              thres = exact_glr_cusum_thres,
                              max_sample = max_sample,
                              num_repeat = num_repeat,
                              run_fn = run_ber_glr_cusum)
WAD_exact_glr_cusum <- WAD_ber_fn(p_pre, p_post,
                              thres = exact_glr_cusum_thres,
                              max_sample = max_sample,
                              num_repeat = num_repeat,
                              run_fn = run_ber_glr_cusum)

# 3. e-SR with unknown p_post
format(Sys.time(), "%c")
set.seed(100)
ARL_eSR_mix <- ARL_ber_fn(p_pre, p_post,
                           thres = log(ARL_target),
                           max_sample = max_sample,
                           num_repeat = num_repeat,
                           run_fn = run_ber_eSR_mix)
WAD_eSR_mix <- WAD_ber_fn(p_pre, p_post,
                           thres = log(ARL_target),
                           max_sample = max_sample,
                           num_repeat = num_repeat,
                           run_fn = run_ber_eSR_mix)

# 4. e-CUSUM with unknown p_post and optimal threshold
format(Sys.time(), "%c")
set.seed(100)
exact_eCS_thres <- find_exact_thres(p_pre, p_post,
                                    ARL_target, max_sample, num_repeat,
                                    run_fn = run_ber_eCS_mix)

ARL_exact_eCS_mix <- ARL_ber_fn(p_pre, p_post,
                          thres = exact_eCS_thres,
                          max_sample = max_sample,
                          num_repeat = num_repeat,
                          run_fn = run_ber_eCS_mix)
WAD_exact_eCS_mix <- WAD_ber_fn(p_pre, p_post,
                          thres = exact_eCS_thres,
                          max_sample = max_sample,
                          num_repeat = num_repeat,
                          run_fn = run_ber_eCS_mix)


# Single run visualization
single_run_plot <- function(x_vec, v = NA, plot_raw = FALSE) {
  exact_cusum <- run_ber_cusum_with_memory(x_vec, p_pre, p_post,
                                           thres = exact_cusum_thres)
  exact_glr_cusum <- run_ber_glr_cusum_with_memory(x_vec, p_pre, p_post,
                                                   thres = exact_glr_cusum_thres)

  eCS_mix <- run_stcp(x_vec, stcp_mix, is_SR_type = FALSE)
  eSR_mix <- run_stcp(x_vec, stcp_mix)

  v_inner <- ifelse(is.na(v), Inf, v)

  if (plot_raw) {
    plot(seq_along(x_vec), exact_cusum$m_vec, type = "l", ylim = c(0, 3 * eSR_mix$stcp_obj$log_one_over_alph), xlab = "n", ylab = "CP-stat")
    abline(v = v_inner, lty = 2)
    abline(h = exact_cusum_thres)
    lines(seq_along(x_vec), exact_glr_cusum$m_vec, col = 2)
    abline(h = exact_glr_cusum_thres, col = 2)
    lines(seq_along(x_vec), eCS_mix$log_mix_e_vec, col = 3)
    abline(h = exact_eCS_thres, col = 3)
    lines(seq_along(x_vec), eSR_mix$log_mix_e_vec, col = 4)
    abline(h = eSR_mix$stcp_obj$log_one_over_alpha, col = 4)
  } else {
    plot(seq_along(x_vec), exact_cusum$m_vec / exact_cusum_thres, type = "l", ylim = c(0, 3), xlab = "n", ylab = "CP-stat")
    abline(v = v_inner, lty = 2)
    abline(h=1)
    lines(seq_along(x_vec), exact_glr_cusum$m_vec / exact_glr_cusum_thres, col = 2)
    lines(seq_along(x_vec), eCS_mix$log_mix_e_vec / exact_eCS_thres, col = 3)
    lines(seq_along(x_vec), eSR_mix$log_mix_e_vec / eSR_mix$stcp_obj$log_one_over_alpha, col = 4)
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
