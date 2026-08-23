################################################################################
# Symbolic Likelihood for Moving-Window Temperature-Humidity Index (THI) Study
# 
# Description: Production-ready replication pipeline for simulation and empirical
#              benchmark experiments in computational statistics (CSDA).
# Output Dir : results_symbolic_thi/ (Data, LaTeX summary tables, Figures)
################################################################################

options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)

# ------------------------------------------------------------------------------
# 0. Global Environment & Configuration
# ------------------------------------------------------------------------------
CONFIG <- list(
  seed                  = 20260818L,
  results_dir           = Sys.getenv("THI_RESULTS_DIR", unset = "results_symbolic_thi"),
  n_cores               = max(1L, parallel::detectCores() - 1L),
  main_reps             = 1000L,
  robustness_reps       = 500L,
  lognormal_reps        = 500L,
  main_models           = c("normal", "logistic", "shifted_gamma"),
  sample_sizes          = c(365L, 1095L, 3650L),
  bins                  = c(10L, 50L, 100L),
  bin_types             = c("EW", "QU"),
  primary_optimizers    = c("BFGS", "L-BFGS-B"),
  robustness_models     = c("normal", "logistic", "shifted_gamma"),
  robustness_bins       = 50L,
  robustness_optimizers = c("Nelder-Mead", "CG"),
  lognormal_model       = "shifted_lognormal",
  empirical_models      = c("normal", "logistic"),
  empirical_optimizers = c("BFGS", "L-BFGS-B"),
  empirical_bins        = 50L,
  window_width          = 3650L,
  step_size             = 182L,
  nasa_start            = as.Date("1984-01-01"),
  nasa_end              = as.Date("2025-12-31"),
  maxit_primary         = 5000L,
  time_K                = 30L,
  time_warm             = 3L
)

# Create execution environment directories
dir.create(CONFIG$results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CONFIG$results_dir, "cache"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Package Installation & Load
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "dplyr", "tibble", "purrr", "readr", "tidyr", "ggplot2",
  "stringr", "jsonlite", "httr2", "future", "future.apply", 
  "parallel", "R.utils", "patchwork"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

set.seed(CONFIG$seed)

# ------------------------------------------------------------------------------
# 2. Labels & Utility Functions
# ------------------------------------------------------------------------------
model_label <- function(model) {
  dplyr::case_when(
    model == "normal"            ~ "Normal",
    model == "logistic"          ~ "Logistic",
    model == "shifted_gamma"     ~ "Shifted Gamma",
    model == "shifted_lognormal" ~ "Shifted Lognormal",
    TRUE ~ model
  )
}

safe_median <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else median(x) }
safe_mean   <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else mean(x) }
safe_ratio  <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_); num / den
}
safe_cor    <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); if (sum(ok) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok]))
}
safe_log    <- function(x) log(pmax(x, .Machine$double.xmin))
make_progress_msg <- function(...) cat(sprintf(...), "\n")

# --- High-Resolution Timing --- #
measure_time <- function(fn, K = CONFIG$time_K, warm = CONFIG$time_warm) {
  if (warm > 0L) for (i in seq_len(warm)) fn()
  st <- base::system.time(for (i in seq_len(K)) fn())
  st[["elapsed"]] / K
}

# ------------------------------------------------------------------------------
# 3. Parametrisation & Truth Functions
# ------------------------------------------------------------------------------
truth_parameters <- function(model) {
  if (model == "normal") return(list(mu = 25, sigma = 5))
  if (model == "logistic") return(list(mu = 25, s = 5 * sqrt(3) / pi))
  if (model == "shifted_gamma") return(list(a = 15, shape = 4, rate = 0.4))
  if (model == "shifted_lognormal") {
    sigma2 <- log(1 + 25/100); sigma <- sqrt(sigma2); mu <- log(10) - sigma2/2
    return(list(a = 15, meanlog = mu, sdlog = sigma))
  }
  stop("Unknown model: ", model)
}

truth_moments <- function(model) {
  tr <- truth_parameters(model)
  if (model == "normal") return(c(mean = tr$mu, sd = tr$sigma))
  if (model == "logistic") return(c(mean = tr$mu, sd = tr$s * pi / sqrt(3)))
  if (model == "shifted_gamma") return(c(mean = tr$a + tr$shape / tr$rate, sd = sqrt(tr$shape) / tr$rate))
  if (model == "shifted_lognormal") {
    em <- exp(tr$meanlog + 0.5*tr$sdlog^2); v <- (exp(tr$sdlog^2)-1) * exp(2*tr$meanlog + tr$sdlog^2)
    return(c(mean = tr$a + em, sd = sqrt(v)))
  }
  stop("Unknown model: ", model)
}

simulate_model <- function(model, n) {
  tr <- truth_parameters(model)
  if (model == "normal") return(rnorm(n, tr$mu, tr$sigma))
  if (model == "logistic") return(rlogis(n, tr$mu, tr$s))
  if (model == "shifted_gamma") return(tr$a + rgamma(n, tr$shape, tr$rate))
  if (model == "shifted_lognormal") return(tr$a + rlnorm(n, tr$meanlog, tr$sdlog))
  stop("Unknown model: ", model)
}

unpack_params <- function(model, eta, a0 = NULL) {
  if (model == "normal") return(list(mu = eta[1], sigma = exp(eta[2])))
  if (model == "logistic") return(list(mu = eta[1], s = exp(eta[2])))
  if (model == "shifted_gamma") {
    if (is.null(a0)) stop("a0 required")
    return(list(shape = exp(eta[1]), rate = exp(eta[2]), a = a0 - exp(eta[3])))
  }
  if (model == "shifted_lognormal") {
    if (is.null(a0)) stop("a0 required")
    return(list(meanlog = eta[1], sdlog = exp(eta[2]), a = a0 - exp(eta[3])))
  }
  stop("Unknown model: ", model)
}

implied_moments_from_pars <- function(model, pars) {
  if (model == "normal") c(mean = pars$mu, sd = pars$sigma)
  else if (model == "logistic") c(mean = pars$mu, sd = pars$s * pi / sqrt(3))
  else if (model == "shifted_gamma") c(mean = pars$a + pars$shape / pars$rate, sd = sqrt(pars$shape) / pars$rate)
  else if (model == "shifted_lognormal") {
    em <- exp(pars$meanlog + 0.5*pars$sdlog^2); v <- (exp(pars$sdlog^2)-1) * exp(2*pars$meanlog + pars$sdlog^2)
    c(mean = pars$a + em, sd = sqrt(v))
  } else stop("Unknown model: ", model)
}

density_model <- function(y, model, pars, log = FALSE) {
  if (model == "normal") return(dnorm(y, pars$mu, pars$sigma, log = log))
  if (model == "logistic") return(dlogis(y, pars$mu, pars$s, log = log))
  if (model == "shifted_gamma") {
    x <- y - pars$a; out <- if (log) rep(-Inf, length(y)) else numeric(length(y))
    ok <- x > 0
    out[ok] <- if (log) dgamma(x[ok], pars$shape, pars$rate, log = TRUE) else dgamma(x[ok], pars$shape, pars$rate)
    return(out)
  }
  if (model == "shifted_lognormal") {
    x <- y - pars$a; out <- if (log) rep(-Inf, length(y)) else numeric(length(y))
    ok <- x > 0
    out[ok] <- if (log) dlnorm(x[ok], pars$meanlog, pars$sdlog, log = TRUE) else dlnorm(x[ok], pars$meanlog, pars$sdlog)
    return(out)
  }
  stop("Unknown model: ", model)
}

cdf_model <- function(y, model, pars) {
  if (model == "normal") return(pnorm(y, pars$mu, pars$sigma))
  if (model == "logistic") return(plogis(y, pars$mu, pars$s))
  if (model == "shifted_gamma") {
    x <- y - pars$a; out <- numeric(length(y)); ok <- x > 0; out[ok] <- pgamma(x[ok], pars$shape, pars$rate)
    return(out)
  }
  if (model == "shifted_lognormal") {
    x <- y - pars$a; out <- numeric(length(y)); ok <- x > 0; out[ok] <- plnorm(x[ok], pars$meanlog, pars$sdlog)
    return(out)
  }
  stop("Unknown model: ", model)
}

# ------------------------------------------------------------------------------
# 4. Initialisation & Parameter Bounds
# ------------------------------------------------------------------------------
initial_eta <- function(model, y, a0 = NULL) {
  m <- mean(y); s <- sd(y); if (!is.finite(s) || s <= 0) s <- 1
  gap <- max(0.10 * s, 1e-3)
  if (model == "normal") return(c(m, log(s)))
  if (model == "logistic") return(c(m, log(max(s * sqrt(3)/pi, 1e-3))))
  if (model %in% c("shifted_gamma","shifted_lognormal")) {
    if (is.null(a0)) stop("a0 required"); shift_guess <- min(a0 - gap, min(y) - 1e-6)
    delta <- max(a0 - shift_guess, 1e-8); z <- y - shift_guess; z[z <= 0] <- min(z[z > 0], na.rm=TRUE)/2
    if (model == "shifted_gamma") {
      mz <- mean(z); vz <- var(z); if (!is.finite(vz) || vz <= 0) vz <- max(1e-3, mz^2/10)
      shape <- max(mz^2/vz, 1e-3); rate <- max(mz/vz, 1e-3); return(c(log(shape), log(rate), log(delta)))
    } else {
      lz <- log(pmax(z, 1e-8)); return(c(mean(lz), log(max(sd(lz), 1e-3)), log(delta)))
    }
  }
  stop("Unknown model: ", model)
}

parameter_bounds <- function(model) {
  min_shape_log <- log(1e-4)
  if (model == "normal") return(list(lower = c(-200, -20), upper = c(200, 20)))
  if (model == "logistic") return(list(lower = c(-200, -20), upper = c(200, 20)))
  if (model == "shifted_gamma") return(list(lower = c(min_shape_log, -10, -20), upper = c(10, 10, 20)))
  if (model == "shifted_lognormal") return(list(lower = c(-20, min_shape_log, -20), upper = c(20, 10, 20)))
  stop("Unknown model: ", model)
}

# ------------------------------------------------------------------------------
# 5. Likelihood Functions & Symbol Generators
# ------------------------------------------------------------------------------
loglik_raw <- function(eta, y, model, a0 = NULL) {
  pars <- unpack_params(model, eta, a0 = a0)
  ll <- sum(density_model(y, model, pars, log = TRUE))
  if (!is.finite(ll)) return(-Inf); ll
}

build_bin_index <- function(y, edges) {
  idx <- findInterval(y, edges, rightmost.closed = TRUE, all.inside = TRUE)
  pmax(1L, pmin(length(edges) - 1L, idx))
}

make_symbol_ew <- function(y, bins, edges = NULL) {
  if (is.null(edges)) {
    rng <- range(y); span <- max(diff(rng), 1e-6); pad <- 0.01 * span + 1e-6
    edges <- seq(rng[1] - pad, rng[2] + pad, length.out = bins + 1L)
  } else bins <- length(edges) - 1L
  idx <- build_bin_index(y, edges); counts <- tabulate(idx, nbins = bins)
  first_occ <- which(counts > 0)[1]
  list(type = "EW", bins = bins, edges = edges, counts = counts, anchor = edges[first_occ])
}

make_symbol_ew_from_counts <- function(counts, edges) {
  first_occ <- which(counts > 0)[1]
  list(type = "EW", bins = length(counts), edges = edges, counts = counts, anchor = edges[first_occ])
}

make_symbol_qu <- function(y, bins) {
  n <- length(y)
  ys <- sort(y)
  probs <- (1:(bins - 1L)) / bins
  indices <- as.integer(round(probs * n))
  indices <- pmin(pmax(indices, 1L), n)
  s <- ys[indices]
  
  if (length(unique(s)) != bins - 1L) {
    s <- s + seq_along(s) * 1e-12
  }
  
  counts <- findInterval(y, c(-Inf, s, Inf))
  counts <- tabulate(counts, nbins = bins)
  
  list(
    type = "QU", bins = bins, internal = s, boundaries = c(-Inf, s, Inf),
    counts = counts, anchor = s[1]
  )
}

loglik_symbolic <- function(eta, symbol, model) {
  pars <- unpack_params(model, eta, a0 = symbol$anchor)
  if (symbol$type == "EW") {
    probs <- cdf_model(symbol$edges[-1], model, pars) - cdf_model(symbol$edges[-length(symbol$edges)], model, pars)
    ll <- sum(symbol$counts * safe_log(probs)); if (!is.finite(ll)) return(-Inf); return(ll)
  }
  if (symbol$type == "QU") {
    probs <- cdf_model(symbol$boundaries[-1], model, pars) - cdf_model(symbol$boundaries[-length(symbol$boundaries)], model, pars)
    dens <- density_model(symbol$internal, model, pars, log = TRUE)
    ll <- sum(dens) + sum(symbol$counts * safe_log(probs)); if (!is.finite(ll)) return(-Inf); return(ll)
  }
  stop("Unknown symbol type")
}

# ------------------------------------------------------------------------------
# 6. Analytical CG Gradients
# ------------------------------------------------------------------------------
grad_ew_normal <- function(eta, sym) {
  mu <- eta[1]; sigma <- exp(eta[2]); F1 <- pnorm(sym$edges, mu, sigma)
  dF <- F1[-1] - F1[-length(F1)]
  mu_term  <- sum(sym$counts * (dnorm(sym$edges[-1], mu, sigma) - dnorm(sym$edges[-length(sym$edges)], mu, sigma)) / pmax(dF, .Machine$double.xmin))
  sig_term <- (sigma) * sum(sym$counts * (((sym$edges[-1]-mu)/sigma) * dnorm(sym$edges[-1], mu, sigma) - ((sym$edges[-length(sym$edges)]-mu)/sigma) * dnorm(sym$edges[-length(sym$edges)], mu, sigma)) / pmax(dF, .Machine$double.xmin))
  c(-mu_term, -sig_term)
}

grad_ew_logistic <- function(eta, sym) {
  mu <- eta[1]; s <- exp(eta[2]); F1 <- plogis(sym$edges, location = mu, scale = s)
  dF <- F1[-1] - F1[-length(F1)]
  pdf <- dlogis(sym$edges, mu, s)
  mu_term  <- sum(sym$counts * (pdf[-1] - pdf[-length(pdf)]) / pmax(dF, .Machine$double.xmin))
  sig_term <- s * sum(sym$counts * (tanh(0.5 * (sym$edges[-1] - mu) / s) * pdf[-1] - tanh(0.5 * (sym$edges[-length(sym$edges)] - mu) / s) * pdf[-length(pdf)]) / pmax(dF, .Machine$double.xmin))
  c(-mu_term, -sig_term)
}

get_analytic_grad <- function(model, sym) {
  if (sym$type != "EW") return(NULL)
  if (model == "normal") return(function(eta) grad_ew_normal(eta, sym))
  if (model == "logistic") return(function(eta) grad_ew_logistic(eta, sym))
  return(NULL)
}

# ------------------------------------------------------------------------------
# 7. Optimisation Wrappers
# ------------------------------------------------------------------------------
optim_control <- function(method) {
  if (method == "L-BFGS-B") list(maxit = CONFIG$maxit_primary, factr = 1e8, trace = 0)
  else if (method == "CG") list(maxit = CONFIG$maxit_primary, reltol = 1e-8, type = 2)
  else list(maxit = CONFIG$maxit_primary, reltol = 1e-8)
}

safe_optim <- function(par, fn, method, lower = NULL, upper = NULL, gr = NULL, retry = TRUE) {
  run_once <- function(start_par) {
    fit <- tryCatch({
      if (method == "L-BFGS-B") {
        optim(par = start_par, fn = fn, method = method, lower = lower, upper = upper, control = optim_control(method))
      } else {
        optim(par = start_par, fn = fn, method = method, gr = gr, control = optim_control(method))
      }
    }, error = function(e) list(par = start_par, value = Inf, convergence = 999, message = conditionMessage(e)))
    if (is.null(fit$message)) fit$message <- NA_character_; fit
  }
  fit1 <- run_once(par)
  ok1 <- is.finite(fit1$value) && fit1$convergence == 0
  if (ok1 || !retry) return(fit1)
  fit2 <- run_once(par + rnorm(length(par), 0, 0.1))
  ok2 <- is.finite(fit2$value) && fit2$convergence == 0
  if (ok2 && (!ok1 || fit2$value < fit1$value)) return(fit2)
  if (is.finite(fit2$value) && fit2$value < fit1$value) return(fit2)
  fit1
}

fit_classical <- function(y, model, method, start_eta = NULL) {
  a0 <- if (model %in% c("shifted_gamma", "shifted_lognormal")) min(y) else NULL
  if (is.null(start_eta)) start_eta <- initial_eta(model, y, a0 = a0)
  bnd <- parameter_bounds(model)
  obj <- function(par) { v <- -loglik_raw(par, y, model, a0 = a0); if (!is.finite(v)) 1e30 else v }
  elapsed <- measure_time(function() safe_optim(start_eta, obj, method, lower = bnd$lower, upper = bnd$upper))
  fit <- safe_optim(start_eta, obj, method, lower = bnd$lower, upper = bnd$upper)
  converged <- is.finite(fit$value) && fit$convergence == 0
  boundary_hit <- any(fit$par <= bnd$lower + 1e-6 | fit$par >= bnd$upper - 1e-6, na.rm = TRUE)
  if (converged) {
    pars <- unpack_params(model, fit$par, a0 = a0)
    mom <- implied_moments_from_pars(model, pars)
    ll <- loglik_raw(fit$par, y, model, a0 = a0)
  } else { mom <- c(mean = NA_real_, sd = NA_real_); ll <- NA_real_ }
  list(converged = converged, eta = fit$par, elapsed = elapsed, value = fit$value,
       loglik = ll, implied = mom, boundary_hit = boundary_hit, message = fit$message)
}

fit_symbolic <- function(y, symbol, model, method, start_eta = NULL) {
  a0 <- symbol$anchor
  if (is.null(start_eta)) start_eta <- initial_eta(model, y, a0 = a0)
  bnd <- parameter_bounds(model)
  obj <- function(par) { v <- -loglik_symbolic(par, symbol, model); if (!is.finite(v)) 1e30 else v }
  gr <- if (method == "CG") get_analytic_grad(model, symbol) else NULL
  elapsed <- measure_time(function() safe_optim(start_eta, obj, method, lower = bnd$lower, upper = bnd$upper, gr = gr))
  fit <- safe_optim(start_eta, obj, method, lower = bnd$lower, upper = bnd$upper, gr = gr)
  converged <- is.finite(fit$value) && fit$convergence == 0
  boundary_hit <- any(fit$par <= bnd$lower + 1e-6 | fit$par >= bnd$upper - 1e-6, na.rm = TRUE)
  if (converged) {
    pars <- unpack_params(model, fit$par, a0 = a0)
    mom <- implied_moments_from_pars(model, pars)
    ll <- loglik_symbolic(fit$par, symbol, model)
  } else { mom <- c(mean = NA_real_, sd = NA_real_); ll <- NA_real_ }
  list(converged = converged, eta = fit$par, elapsed = elapsed, value = fit$value,
       loglik = ll, implied = mom, boundary_hit = boundary_hit, message = fit$message)
}

# ------------------------------------------------------------------------------
# 8. Parallel Simulation Engine
# ------------------------------------------------------------------------------
run_simulation_experiment <- function(models, ns, bins_vec, bin_types, optimizers, reps, study_name) {
  out <- list(); kk <- 1L
  workers <- min(CONFIG$n_cores, reps)
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  for (model in models) for (n in ns) {
    trm <- truth_moments(model)
    make_progress_msg("[%s] model=%s, n=%d | Running %d reps across %d parallel cores...",
                      study_name, model, n, reps, workers)
    rep_results <- future.apply::future_lapply(seq_len(reps), function(rep_id) {
      y <- simulate_model(model, n)
      classical_fits <- setNames(vector("list", length(optimizers)), optimizers)
      for (opt in optimizers) classical_fits[[opt]] <- fit_classical(y, model, opt)
      symbols <- list()
      for (b in bins_vec) {
        if ("EW" %in% bin_types) symbols[[paste0("EW_", b)]] <- make_symbol_ew(y, bins = b)
        if ("QU" %in% bin_types) symbols[[paste0("QU_", b)]] <- make_symbol_qu(y, bins = b)
      }
      sub_out <- list(); sub_kk <- 1L
      for (b in bins_vec) for (bt in bin_types) {
        sym <- symbols[[paste0(bt, "_", b)]]
        for (opt in optimizers) {
          cfit <- classical_fits[[opt]]; sfit <- fit_symbolic(y, sym, model, opt)
          sub_out[[sub_kk]] <- tibble::tibble(
            study = study_name, replicate = rep_id, model = model, optimizer = opt,
            n = n, bins = b, bin_type = bt,
            truth_mean = trm["mean"], truth_sd = trm["sd"],
            class_mean = cfit$implied["mean"], class_sd = cfit$implied["sd"],
            sym_mean = sfit$implied["mean"], sym_sd = sfit$implied["sd"],
            time_classical = cfit$elapsed, time_symbolic = sfit$elapsed,
            class_loglik = cfit$loglik, sym_loglik = sfit$loglik,
            class_conv = cfit$converged, sym_conv = sfit$converged,
            paired_conv = cfit$converged & sfit$converged,
            class_boundary_hit = cfit$boundary_hit, sym_boundary_hit = sfit$boundary_hit
          ); sub_kk <- sub_kk + 1L
        }
      }
      dplyr::bind_rows(sub_out)
    }, future.seed = TRUE)
    out[[kk]] <- dplyr::bind_rows(rep_results); kk <- kk + 1L
  }
  dplyr::bind_rows(out)
}

summarise_simulation <- function(df) {
  df %>% mutate(class_sq_err = (class_mean - truth_mean)^2 + (class_sd - truth_sd)^2,
                sym_sq_err   = (sym_mean   - truth_mean)^2 + (sym_sd   - truth_sd)^2) %>%
    group_by(model, optimizer, n, bins, bin_type) %>% summarise(
      reps = n(), usable = sum(paired_conv, na.rm = TRUE),
      convergence_rate = mean(paired_conv, na.rm = TRUE),
      classical_conv_rate = mean(class_conv, na.rm = TRUE),
      symbolic_conv_rate = mean(sym_conv, na.rm = TRUE),
      median_time_classical = safe_median(time_classical[paired_conv]),
      median_time_symbolic  = safe_median(time_symbolic[paired_conv]),
      speedup = safe_ratio(median_time_classical, median_time_symbolic),
      mse_class = safe_mean(class_sq_err[paired_conv]),
      mse_symbolic = safe_mean(sym_sq_err[paired_conv]),
      rmse_class = sqrt(mse_class), rmse_symbolic = sqrt(mse_symbolic),
      relative_rmse = ifelse(is.finite(mse_class) & mse_class > 0,
                             sqrt(mse_symbolic / mse_class), NA_real_),
      rmse_mean_ratio = safe_ratio(sqrt(safe_mean((sym_mean[paired_conv] - truth_mean[paired_conv])^2)),
                                   sqrt(safe_mean((class_mean[paired_conv] - truth_mean[paired_conv])^2))),
      rmse_sd_ratio   = safe_ratio(sqrt(safe_mean((sym_sd[paired_conv] - truth_sd[paired_conv])^2)),
                                   sqrt(safe_mean((class_sd[paired_conv] - truth_sd[paired_conv])^2))),
      boundary_hit_rate = mean(class_boundary_hit | sym_boundary_hit, na.rm = TRUE)
    ) %>% ungroup() %>% arrange(model, optimizer, n, bins, bin_type)
}

write_sim_latex_rows <- function(df, file) {
  lines <- df %>% mutate(model = model_label(model),
                         line = sprintf("%s & %s & %d & %d & %s & %.4f & %.2f & %.3f \\\\",
                                        model, optimizer, n, bins, bin_type, relative_rmse, speedup, convergence_rate)) %>% pull(line)
  writeLines(lines, con = file)
}

# ------------------------------------------------------------------------------
# 9. NASA POWER API Importer & Empirical Preparation
# ------------------------------------------------------------------------------
fetch_nasa_power_daily <- function(location, latitude, longitude, start_date, end_date, cache_dir) {
  cache_file <- file.path(cache_dir, paste0(str_replace_all(tolower(location),"[^a-z0-9]+","_"), "_nasa_power.csv"))
  if (file.exists(cache_file)) return(readr::read_csv(cache_file, show_col_types = FALSE))
  req <- httr2::request("https://power.larc.nasa.gov/api/temporal/daily/point") |>
    httr2::req_url_query(parameters = "T2M,RH2M", community = "RE",
                         longitude = sprintf("%.4f", longitude),
                         latitude = sprintf("%.4f", latitude),
                         start = format(start_date, "%Y%m%d"),
                         end = format(end_date, "%Y%m%d"),
                         format = "JSON", `time-standard` = "UTC")
  resp <- NULL; max_retries <- 5L
  for (attempt in seq_len(max_retries)) {
    resp <- tryCatch(httr2::req_perform(req),
                     error = function(e) { 
                       make_progress_msg("NASA POWER attempt %d/%d failed for %s: %s", attempt, max_retries, location, e$message)
                       if (attempt < max_retries) Sys.sleep(2*attempt); NULL 
                     })
    if (!is.null(resp)) break
  }
  if (is.null(resp)) stop("Failed to fetch NASA POWER data for ", location)
  js <- httr2::resp_body_json(resp, simplifyVector = TRUE)
  df <- tibble(date = as.Date(names(js$properties$parameter$T2M), "%Y%m%d"),
               T2M  = as.numeric(js$properties$parameter$T2M),
               RH2M = as.numeric(js$properties$parameter$RH2M)) |>
    filter(is.finite(T2M), is.finite(RH2M), T2M > -900, RH2M > -900) |> arrange(date)
  readr::write_csv(df, cache_file); Sys.sleep(0.5); df
}

compute_thi <- function(T, RH) T - (0.55 - 0.0055 * RH) * (T - 14.5)

location_table <- tibble::tribble(
  ~location, ~latitude, ~longitude, ~climate_classification,
  "Yakutsk",      62.0355,  129.6755, "Subarctic",
  "Singapore",     1.3521,  103.8198, "Tropical rainforest",
  "Death Valley", 36.4600, -116.8650, "Desert",
  "Melbourne",   -37.8136,  144.9631, "Oceanic",
  "Chicago",      41.8781,  -87.6298, "Humid continental"
)

# ------------------------------------------------------------------------------
# 10. Empirical Moving-Window Engine
# ------------------------------------------------------------------------------
run_empirical_one <- function(df_loc, location_name, model, optimizer, bins, window_width, step_size) {
  y <- df_loc$THI; dates <- df_loc$date; n_total <- length(y)
  if (n_total < window_width) stop("Not enough observations for ", location_name)
  rng <- range(y, na.rm = TRUE); span <- max(diff(rng), 1e-6); pad <- 0.01 * span + 1e-6
  fixed_edges <- seq(rng[1] - pad, rng[2] + pad, length.out = bins + 1L); nb <- length(fixed_edges) - 1L
  starts <- seq.int(1L, n_total - window_width + 1L, by = step_size)
  bin_ids <- build_bin_index(y, fixed_edges)
  s <- starts[1]; current_counts <- tabulate(bin_ids[s:(s + window_width - 1L)], nbins = nb)
  prev_class_eta <- prev_sym_eta <- NULL
  out <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    s <- starts[i]; e <- s + window_width - 1L
    if (i > 1) {
      prev_s <- starts[i-1L]; prev_e <- prev_s + window_width - 1L
      current_counts <- current_counts -
        tabulate(bin_ids[prev_s:(s-1L)], nbins = nb) +
        tabulate(bin_ids[(prev_e+1L):e],  nbins = nb)
    }
    y_win <- y[s:e]; symbol <- make_symbol_ew_from_counts(current_counts, fixed_edges)
    cfit <- fit_classical(y_win, model, optimizer, start_eta = prev_class_eta)
    sfit <- fit_symbolic(y_win, symbol, model, optimizer, start_eta = prev_sym_eta)
    prev_class_eta <- if (cfit$converged) cfit$eta else initial_eta(model, y_win, a0 = if (grepl("shifted", model)) min(y_win) else NULL)
    prev_sym_eta   <- if (sfit$converged) sfit$eta else initial_eta(model, y_win, a0 = symbol$anchor)
    out[[i]] <- tibble(location = location_name, model = model, optimizer = optimizer,
                       window_index = i, start_date = dates[s], end_date = dates[e],
                       class_mean = cfit$implied["mean"], class_sd = cfit$implied["sd"],
                       sym_mean = sfit$implied["mean"], sym_sd = sfit$implied["sd"],
                       time_classical = cfit$elapsed, time_symbolic = sfit$elapsed,
                       class_conv = cfit$converged, sym_conv = sfit$converged,
                       paired_conv = cfit$converged & sfit$converged)
  }
  bind_rows(out)
}

run_empirical_study <- function() {
  cache_dir <- file.path(CONFIG$results_dir, "cache")
  detail_all <- list(); kk <- 1L
  for (i in seq_len(nrow(location_table))) {
    loc <- location_table[i, ]
    make_progress_msg("[EMPIRICAL] Downloading/reading %s ...", loc$location)
    df_raw <- fetch_nasa_power_daily(loc$location, loc$latitude, loc$longitude,
                                     CONFIG$nasa_start, CONFIG$nasa_end, cache_dir)
    df_loc <- df_raw %>% mutate(THI = compute_thi(T2M, RH2M)) %>% filter(is.finite(THI)) %>% arrange(date)
    readr::write_csv(df_loc, file.path(CONFIG$results_dir,
                                       paste0(str_replace_all(tolower(loc$location), "[^a-z0-9]+","_"),"_thi_series.csv")))
    for (model in CONFIG$empirical_models) for (opt in CONFIG$empirical_optimizers) {
      make_progress_msg("[EMPIRICAL] %s | %s | %s", loc$location, model, opt)
      detail_all[[kk]] <- run_empirical_one(df_loc, loc$location, model, opt,
                                            CONFIG$empirical_bins, CONFIG$window_width, CONFIG$step_size)
      kk <- kk + 1L
    }
  }
  bind_rows(detail_all)
}

summarise_empirical <- function(df) {
  df %>% group_by(location, model, optimizer) %>% summarise(
    N_windows = n(), usable = sum(paired_conv, na.rm = TRUE),
    median_time_classical = safe_median(time_classical[paired_conv]),
    median_time_symbolic  = safe_median(time_symbolic[paired_conv]),
    speedup = safe_ratio(median_time_classical, median_time_symbolic),
    correlation = safe_cor(class_mean[paired_conv], sym_mean[paired_conv]),
    sd_correlation = safe_cor(class_sd[paired_conv], sym_sd[paired_conv]),
    rmse_discrepancy = sqrt(safe_mean((sym_mean[paired_conv]-class_mean[paired_conv])^2 +
                                        (sym_sd[paired_conv]  -class_sd[paired_conv])^2)),
    convergence_rate = mean(paired_conv, na.rm = TRUE),
    classical_conv_rate = mean(class_conv, na.rm = TRUE),
    symbolic_conv_rate = mean(sym_conv, na.rm = TRUE)
  ) %>% ungroup() %>% arrange(location, model, optimizer)
}

write_empirical_latex_rows <- function(df, file) {
  lines <- df %>% mutate(model = model_label(model),
                         line = sprintf("%s & %s & %s & %d & %.4f & %.4f & %.2f & %.3f \\\\",
                                        location, model, optimizer, N_windows,
                                        median_time_classical, median_time_symbolic,
                                        speedup, correlation)) %>% pull(line)
  writeLines(lines, con = file)
}

# ------------------------------------------------------------------------------
# 11. Plotting Generators
# ------------------------------------------------------------------------------
plot_main_speedup <- function(df, file) {
  p <- df %>% mutate(model = factor(model_label(model), levels = c("Normal","Logistic","Shifted Gamma")),
                     n = factor(n), bin_type = factor(bin_type, levels = c("EW","QU"))) %>%
    ggplot(aes(x = n, y = speedup, colour = factor(bins), shape = bin_type,
               group = interaction(bins, bin_type))) +
    geom_line() + geom_point(size = 2.2) + facet_grid(model ~ optimizer) +
    labs(x = "Sample size n", y = "Speed-up (classical / symbolic)",
         colour = "B-1", shape = "Bin type") + theme_bw(base_size = 12)
  ggsave(file, p, width = 13, height = 7.5, dpi = 300)
}

plot_main_relative_rmse <- function(df, file) {
  p <- df %>% mutate(model = factor(model_label(model), levels = c("Normal","Logistic","Shifted Gamma")),
                     n = factor(n), bin_type = factor(bin_type, levels = c("EW","QU"))) %>%
    ggplot(aes(x = n, y = relative_rmse, colour = factor(bins), shape = bin_type,
               group = interaction(bins, bin_type))) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_line() + geom_point(size = 2.2) + facet_grid(model ~ optimizer) +
    labs(x = "Sample size n", y = "Relative RMSE (symbolic / classical)",
         colour = "B-1", shape = "Bin type") + theme_bw(base_size = 12)
  ggsave(file, p, width = 13, height = 7.5, dpi = 300)
}

plot_empirical_overlay <- function(file_path = file.path(CONFIG$results_dir, "fig_empirical_thi_bar_overlay.png")) {
  bins <- CONFIG$empirical_bins
  loc_files <- list.files(CONFIG$results_dir, pattern = "_thi_series\\.csv$", full.names = TRUE)
  if (length(loc_files) == 0) return(NULL)
  
  emp_data <- purrr::map_dfr(loc_files, function(f) {
    loc_key <- stringr::str_remove(basename(f), "_thi_series\\.csv$")
    match_idx <- match(loc_key, stringr::str_replace_all(tolower(location_table$location), "[^a-z0-9]+", "_"))
    loc_clean <- if (!is.na(match_idx)) location_table$location[match_idx] else stringr::str_to_title(loc_key)
    readr::read_csv(f, show_col_types = FALSE) %>% dplyr::mutate(location = loc_clean)
  })
  
  emp_data <- emp_data %>%
    dplyr::left_join(location_table %>% dplyr::select(location, climate_classification), by = "location") %>%
    dplyr::mutate(
      facet_label = paste0(location, " (", climate_classification, ")"),
      location = factor(location, levels = location_table$location)
    )
  
  facet_levels <- unique(emp_data$facet_label[order(match(emp_data$location, location_table$location))])
  emp_data$facet_label <- factor(emp_data$facet_label, levels = facet_levels)
  
  snapshot_windows <- tibble::tribble(
    ~period_label, ~start_yr, ~end_yr,
    "1985–1995",   1985,      1995,
    "1995–2005",   1995,      2005,
    "2005–2015",   2005,      2015,
    "2015–2025",   2015,      2025
  )
  
  hist_bars_list <- list()
  for (loc_name in levels(emp_data$location)) {
    sub_df <- emp_data %>% dplyr::filter(location == loc_name)
    facet_lbl <- sub_df$facet_label[1]
    y_all <- sub_df$THI
    
    rng <- range(y_all, na.rm = TRUE)
    edges <- seq(rng[1], rng[2], length.out = bins + 1L)
    bin_width <- diff(edges)[1]
    
    for (w in seq_len(nrow(snapshot_windows))) {
      win_info <- snapshot_windows[w, ]
      win_df <- sub_df %>% 
        dplyr::filter(date >= as.Date(paste0(win_info$start_yr, "-01-01")) & 
                        date <  as.Date(paste0(win_info$end_yr, "-01-01")))
      
      if (nrow(win_df) == 0) next
      
      counts <- tabulate(build_bin_index(win_df$THI, edges), nbins = bins)
      max_c <- max(counts, na.rm = TRUE)
      center_date <- as.Date(paste0(as.integer((win_info$start_yr + win_info$end_yr)/2), "-07-01"))
      
      for (b in seq_along(counts)) {
        if (counts[b] == 0) next
        bar_height <- bin_width * 0.90
        y_bottom <- edges[b]
        y_top <- y_bottom + bar_height
        time_span_days <- (counts[b] / max_c) * (365 * 4.5)
        
        hist_bars_list[[length(hist_bars_list) + 1L]] <- tibble::tibble(
          location = loc_name,
          facet_label = facet_lbl,
          period = win_info$period_label,
          xmin = center_date,
          xmax = center_date + time_span_days,
          ymin = y_bottom,
          ymax = y_top
        )
      }
    }
  }
  
  bars_df <- dplyr::bind_rows(hist_bars_list)
  bars_df$facet_label <- factor(bars_df$facet_label, levels = facet_levels)
  
  period_colors <- c(
    "1985–1995" = "#059669",
    "1995–2005" = "#2563EB",
    "2005–2015" = "#D97706",
    "2015–2025" = "#DC2626"
  )
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = emp_data, 
      ggplot2::aes(x = date, y = THI), 
      color = "#334155", linewidth = 0.22, alpha = 0.40
    ) +
    ggplot2::geom_rect(
      data = bars_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = period),
      color = "#FFFFFF", linewidth = 0.12, alpha = 0.75
    ) +
    ggplot2::facet_wrap(~facet_label, scales = "free_y", ncol = 1) +
    ggplot2::scale_fill_manual(
      values = period_colors, 
      name = "10-Year Decadal Baseline Period"
    ) +
    ggplot2::scale_x_date(
      date_breaks = "5 years", 
      date_labels = "%Y", 
      expand = c(0.01, 0)
    ) +
    ggplot2::labs(
      x = "Year", 
      y = "Daily Temperature-Humidity Index (THI)"
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 10, r = 12, b = 10, l = 12),
      axis.title.x = ggplot2::element_text(face = "bold", size = 14, color = "#1E293B", margin = ggplot2::margin(t = 8, b = 8)),
      axis.title.y = ggplot2::element_text(face = "bold", size = 14, color = "#1E293B", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 11, color = "#1E293B", face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", size = 13, color = "#0F172A", hjust = 0.02),
      strip.background = ggplot2::element_rect(fill = "#F1F5F9", color = "#CBD5E1", linewidth = 0.6),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#E2E8F0", linewidth = 0.35),
      panel.border = ggplot2::element_rect(color = "#94A3B8", fill = NA, linewidth = 0.6),
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = ggplot2::element_text(face = "bold", size = 12, color = "#0F172A"),
      legend.text = ggplot2::element_text(size = 11, face = "bold", color = "#334155"),
      legend.key.size = grid::unit(0.8, "cm"),
      legend.margin = ggplot2::margin(t = 4, b = 2)
    )
  
  ggplot2::ggsave(file_path, plot = p, width = 15, height = 13, dpi = 300)
  make_progress_msg("Saved publication-grade discrete histogram overlay plot to: %s", file_path)
  return(p)
}

# ------------------------------------------------------------------------------
# 12. Multi-Start Diagnostic
# ------------------------------------------------------------------------------
run_multistart_check <- function() {
  cat("\nRunning multi-start diagnostic: Shifted Gamma, n=3650, B-1=50...\n")
  set.seed(CONFIG$seed)
  y <- simulate_model("shifted_gamma", 3650L); sym <- make_symbol_ew(y, bins = 50L)
  base_eta <- initial_eta("shifted_gamma", y, a0 = sym$anchor)
  base_fit <- fit_symbolic(y, sym, "shifted_gamma", "L-BFGS-B", start_eta = base_eta)
  pars_def <- unpack_params("shifted_gamma", base_fit$eta, a0 = sym$anchor)
  mean_def <- pars_def$a + pars_def$shape / pars_def$rate
  successes <- 0L
  for (i in 1:20) {
    s_eta <- base_eta + rnorm(length(base_eta), 0, 0.05)
    fit <- fit_symbolic(y, sym, "shifted_gamma", "L-BFGS-B", start_eta = s_eta)
    if (fit$converged && is.finite(fit$value)) {
      p <- unpack_params("shifted_gamma", fit$eta, a0 = sym$anchor)
      if (abs((p$a + p$shape/p$rate) - mean_def) < 0.1) successes <- successes + 1L
    }
  }
  cat(sprintf("Multi-start recovery rate: %d/20 (%.1f%%)\n", successes, 100*successes/20))
  successes / 20
}

# ------------------------------------------------------------------------------
# 13. Master Pipeline Orchestration
# ------------------------------------------------------------------------------
run_all_analysis <- function() {
  make_progress_msg("====================================================")
  make_progress_msg("Symbolic THI Analysis Pipeline")
  make_progress_msg("Results directory: %s", CONFIG$results_dir)
  make_progress_msg("Parallel cores: %d", CONFIG$n_cores)
  make_progress_msg("====================================================")
  
  writeLines(jsonlite::toJSON(CONFIG, auto_unbox = TRUE, pretty = TRUE),
             file.path(CONFIG$results_dir, "run_config.json"))
  writeLines(capture.output(sessionInfo()),
             file.path(CONFIG$results_dir, "session_info.txt"))
  
  # 1. Primary Simulation Study
  main_detail <- run_simulation_experiment(CONFIG$main_models, CONFIG$sample_sizes,
                                           CONFIG$bins, CONFIG$bin_types,
                                           CONFIG$primary_optimizers, CONFIG$main_reps,
                                           "main_primary")
  main_summary <- summarise_simulation(main_detail)
  readr::write_csv(main_detail, file.path(CONFIG$results_dir, "sim_main_primary_detailed.csv"))
  readr::write_csv(main_summary, file.path(CONFIG$results_dir, "sim_main_primary_summary.csv"))
  write_sim_latex_rows(main_summary, file.path(CONFIG$results_dir, "sim_main_primary_summary_rows.tex"))
  plot_main_speedup(main_summary,           file.path(CONFIG$results_dir, "fig_speedup.png"))
  plot_main_relative_rmse(main_summary,     file.path(CONFIG$results_dir, "fig_relrmse.png"))
  
  # 2. Optimiser Robustness Study
  robustness_detail <- run_simulation_experiment(CONFIG$robustness_models, CONFIG$sample_sizes,
                                                 CONFIG$robustness_bins, CONFIG$bin_types,
                                                 CONFIG$robustness_optimizers, CONFIG$robustness_reps,
                                                 "optimizer_robustness")
  robustness_summary <- summarise_simulation(robustness_detail)
  readr::write_csv(robustness_detail, file.path(CONFIG$results_dir, "sim_optimizer_robustness_detailed.csv"))
  readr::write_csv(robustness_summary, file.path(CONFIG$results_dir, "sim_optimizer_robustness_summary.csv"))
  write_sim_latex_rows(robustness_summary, file.path(CONFIG$results_dir, "sim_optimizer_robustness_summary_rows.tex"))
  
  # 3. Shifted-Lognormal Sensitivity Analysis
  lognormal_detail <- run_simulation_experiment(CONFIG$lognormal_model, CONFIG$sample_sizes,
                                                CONFIG$bins, CONFIG$bin_types,
                                                CONFIG$primary_optimizers, CONFIG$lognormal_reps,
                                                "shifted_lognormal_sensitivity")
  lognormal_summary <- summarise_simulation(lognormal_detail)
  readr::write_csv(lognormal_detail, file.path(CONFIG$results_dir, "sim_shifted_lognormal_detailed.csv"))
  readr::write_csv(lognormal_summary, file.path(CONFIG$results_dir, "sim_shifted_lognormal_summary.csv"))
  write_sim_latex_rows(lognormal_summary, file.path(CONFIG$results_dir, "sim_shifted_lognormal_summary_rows.tex"))
  
  # 4. Empirical Moving-Window Application
  empirical_detail <- run_empirical_study()
  empirical_summary <- summarise_empirical(empirical_detail)
  readr::write_csv(empirical_detail, file.path(CONFIG$results_dir, "empirical_detailed.csv"))
  readr::write_csv(empirical_summary, file.path(CONFIG$results_dir, "empirical_summary.csv"))
  write_empirical_latex_rows(empirical_summary, file.path(CONFIG$results_dir, "empirical_summary_rows.tex"))
  
  # 5. Publication Overlay Figure Generation
  plot_empirical_overlay()
  
  # 6. Diagnostic Check
  run_multistart_check()
  
  cat("\n--- MAIN SIMULATION SUMMARY ---\n");         print(main_summary, n = Inf)
  cat("\n--- OPTIMISER ROBUSTNESS SUMMARY ---\n");    print(robustness_summary, n = Inf)
  cat("\n--- SHIFTED-LOGNORMAL SENSITIVITY ---\n");   print(lognormal_summary, n = Inf)
  cat("\n--- EMPIRICAL MOVING-WINDOW SUMMARY ---\n"); print(empirical_summary, n = Inf)
  cat("\nAll output files saved to: ", normalizePath(CONFIG$results_dir, winslash="/"), "\n")
  
  invisible(list(main_detail=main_detail, main_summary=main_summary,
                 robustness_detail=robustness_detail, robustness_summary=robustness_summary,
                 lognormal_detail=lognormal_detail, lognormal_summary=lognormal_summary,
                 empirical_detail=empirical_detail, empirical_summary=empirical_summary))
}

# ------------------------------------------------------------------------------
# 14. Entry Point / Execution
# ------------------------------------------------------------------------------
results <- run_all_analysis()