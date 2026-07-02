################################################################################
#####                                Utils                                 #####
################################################################################

loss_mse <- function(h_proxy, h_fore)  (h_proxy - h_fore)^2
loss_qlike <- function(h_proxy, h_fore) h_proxy/h_fore - log(h_proxy/h_fore) - 1
loss_mse_log <- function(h_proxy, h_fore) (log(h_proxy) - log(h_fore))^2
loss_mse_sd <- function(h_proxy, h_fore) (sqrt(h_proxy) - sqrt(h_fore))^2
loss_mse_prop <- function(h_proxy, h_fore) (h_proxy / h_fore - 1)^2
loss_mae <- function(h_proxy, h_fore) abs(h_proxy - h_fore)
loss_mae_log <- function(h_proxy, h_fore) abs(log(h_proxy) - log(h_fore))
loss_mae_sd <- function(h_proxy, h_fore) abs(sqrt(h_proxy) - sqrt(h_fore))
loss_mae_prop <- function(h_proxy, h_fore) abs(h_proxy / h_fore - 1)
compute_bias <- function(h_proxy, h_fore) mean(h_fore - h_proxy, na.rm = TRUE)

table_bias <- function(file) {
  out_500 <- rbind(
    compute_bias(file$sigma, file$garch.n.500),
    compute_bias(file$sigma, file$garch.t.500),
    compute_bias(file$sigma, file$figarch.n.500),
    compute_bias(file$sigma, file$figarch.t.500),
    compute_bias(file$sigma, file$gas.n.500),
    compute_bias(file$sigma, file$gas.t.500),
    compute_bias(file$sigma, file$svb.n.500),
    compute_bias(file$sigma, file$svb.t.500),
    compute_bias(file$sigma, file$ms.n.500),
    compute_bias(file$sigma, file$ms.t.500)
    )
  out_1000 <- rbind(
    compute_bias(file$sigma, file$garch.n.1000),
    compute_bias(file$sigma, file$garch.t.1000),
    compute_bias(file$sigma, file$figarch.n.1000),
    compute_bias(file$sigma, file$figarch.t.1000),
    compute_bias(file$sigma, file$gas.n.1000),
    compute_bias(file$sigma, file$gas.t.1000),
    compute_bias(file$sigma, file$svb.n.1000),
    compute_bias(file$sigma, file$svb.t.1000),
    compute_bias(file$sigma, file$ms.n.1000),
    compute_bias(file$sigma, file$ms.t.1000)
  )
  out_2500 <- rbind(
    compute_bias(file$sigma, file$garch.n.2500),
    compute_bias(file$sigma, file$garch.t.2500),
    compute_bias(file$sigma, file$figarch.n.2500),
    compute_bias(file$sigma, file$figarch.t.2500),
    compute_bias(file$sigma, file$gas.n.2500),
    compute_bias(file$sigma, file$gas.t.2500),
    compute_bias(file$sigma, file$svb.n.2500),
    compute_bias(file$sigma, file$svb.t.2500),
    compute_bias(file$sigma, file$ms.n.2500),
    compute_bias(file$sigma, file$ms.t.2500)
  )
  out <- round(cbind(out_500, out_1000, out_2500), 4)
  return(out)
}

table_loss <- function(file, w = 0.00, loss = "MSE") {
  
  metric <- switch(loss,
    "MSE" = loss_mse,
    "QLIKE" = loss_qlike,
    "MSE-LOG" = loss_mse_log,
    "MSE-SD" = loss_mse_sd ,
    "MSE-PROP" = loss_mse_prop,
    "MAE" = loss_mae,
    "MAE-LOG" = loss_mae_log,
    "MAE-SD" = loss_mae_sd,
    "MAE-PROP" = loss_mae_prop,
    stop("Error."))
  
  pre_out_500 <- cbind(
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.t.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.t.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.t.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.t.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.t.500))
  
  pre_out_1000 <- cbind(
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.t.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.t.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.t.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.t.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.t.1000))
  
  pre_out_2500 <- cbind(
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.t.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.t.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.t.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.t.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.t.2500))
  
  out_500  <- apply(pre_out_500,  2, mean)
  out_1000 <- apply(pre_out_1000, 2, mean)
  out_2500 <- apply(pre_out_2500, 2, mean)
  
  out <- cbind(out_500, out_1000, out_2500)
  return(out)
}

table_mcs_loss <- function(file, w = 0.00, loss = "MSE") {
  
  metric <- switch(loss,
    "MSE" = loss_mse,
    "QLIKE" = loss_qlike,
    "MSE-LOG" = loss_mse_log,
    "MSE-SD" = loss_mse_sd ,
    "MSE-PROP" = loss_mse_prop,
    "MAE" = loss_mae,
    "MAE-LOG" = loss_mae_log,
    "MAE-SD" = loss_mae_sd,
    "MAE-PROP" = loss_mae_prop,
    stop("Error."))
  
  pre_out_500 <- cbind(
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.t.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.t.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.t.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.t.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.n.500),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.t.500))
  
  pre_out_1000 <- cbind(
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.t.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.t.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.t.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.t.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.n.1000),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.t.1000))
  
  pre_out_2500 <- cbind(
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$garch.t.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$figarch.t.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$gas.t.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$svb.t.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.n.2500),
    metric(file$sigma * (1 -w) + file$returns * w, file$ms.t.2500))
  
  mcs_0500 <- mcs_1000 <- mcs_2500 <- rep(0, 10)
  
  mcs_0500[estMCS.quick(pre_out_500 , test = "t.range", B = 5000, l = 1, alpha = 0.25)] <- 1
  mcs_1000[estMCS.quick(pre_out_1000, test = "t.range", B = 5000, l = 1, alpha = 0.25)] <- 1
  mcs_2500[estMCS.quick(pre_out_2500, test = "t.range", B = 5000, l = 1, alpha = 0.25)] <- 1

  out <- cbind(mcs_0500, mcs_1000, mcs_2500)
  return(out)
}
  
table_latex <- function(loss_func, w) {
  loss_dgp_garch   <- cbind(table_loss(garch_n, w, loss = loss_func),   table_loss(garch_n_out_m, w, loss = loss_func),   table_loss(garch_n_out_w, w, loss = loss_func),   table_loss(garch_t, w, loss = loss_func),   table_loss(garch_t_out_m, w, loss = loss_func),   table_loss(garch_t_out_w, w, loss = loss_func))
  loss_dgp_figarch <- cbind(table_loss(figarch_n, w, loss = loss_func), table_loss(figarch_n_out_m, w, loss = loss_func), table_loss(figarch_n_out_w, w, loss = loss_func), table_loss(figarch_t, w, loss = loss_func), table_loss(figarch_t_out_m, w, loss = loss_func), table_loss(figarch_t_out_w, w, loss = loss_func))
  loss_dgp_dcs     <- cbind(table_loss(dcs_n, w, loss = loss_func),     table_loss(dcs_n_out_m, w, loss = loss_func),     table_loss(dcs_n_out_w, w, loss = loss_func),     table_loss(dcs_t, w, loss = loss_func),     table_loss(dcs_t_out_m, w, loss = loss_func),     table_loss(dcs_t_out_w, w, loss = loss_func))
  loss_dgp_sv      <- cbind(table_loss(sv_n, w, loss = loss_func),      table_loss(sv_n_out_m, w, loss = loss_func),      table_loss(sv_n_out_w, w, loss = loss_func),      table_loss(sv_t, w, loss = loss_func),      table_loss(sv_t_out_m, w, loss = loss_func),      table_loss(sv_t_out_w, w, loss = loss_func))
  loss_dgp_ms      <- cbind(table_loss(ms_n, w, loss = loss_func),      table_loss(ms_n_out_m, w, loss = loss_func),      table_loss(ms_n_out_w, w, loss = loss_func),      table_loss(ms_t, w, loss = loss_func),      table_loss(ms_t_out_m, w, loss = loss_func),      table_loss(ms_t_out_w, w, loss = loss_func))
  
  loss_mcs_dgp_garch   <- cbind(table_mcs_loss(garch_n, w, loss = loss_func),   table_mcs_loss(garch_n_out_m, w, loss = loss_func),   table_mcs_loss(garch_n_out_w, w, loss = loss_func),   table_mcs_loss(garch_t, w, loss = loss_func),   table_mcs_loss(garch_t_out_m, w, loss = loss_func),   table_mcs_loss(garch_t_out_w, w, loss = loss_func))
  loss_mcs_dgp_figarch <- cbind(table_mcs_loss(figarch_n, w, loss = loss_func), table_mcs_loss(figarch_n_out_m, w, loss = loss_func), table_mcs_loss(figarch_n_out_w, w, loss = loss_func), table_mcs_loss(figarch_t, w, loss = loss_func), table_mcs_loss(figarch_t_out_m, w, loss = loss_func), table_mcs_loss(figarch_t_out_w, w, loss = loss_func))
  loss_mcs_dgp_dcs     <- cbind(table_mcs_loss(dcs_n, w, loss = loss_func),     table_mcs_loss(dcs_n_out_m, w, loss = loss_func),     table_mcs_loss(dcs_n_out_w, w, loss = loss_func),     table_mcs_loss(dcs_t, w, loss = loss_func),     table_mcs_loss(dcs_t_out_m, w, loss = loss_func),     table_mcs_loss(dcs_t_out_w, w, loss = loss_func))
  loss_mcs_dgp_sv      <- cbind(table_mcs_loss(sv_n, w, loss = loss_func),      table_mcs_loss(sv_n_out_m, w, loss = loss_func),      table_mcs_loss(sv_n_out_w, w, loss = loss_func),      table_mcs_loss(sv_t, w, loss = loss_func),      table_mcs_loss(sv_t_out_m, w, loss = loss_func),      table_mcs_loss(sv_t_out_w, w, loss = loss_func))
  loss_mcs_dgp_ms      <- cbind(table_mcs_loss(ms_n, w, loss = loss_func),      table_mcs_loss(ms_n_out_m, w, loss = loss_func),      table_mcs_loss(ms_n_out_w, w, loss = loss_func),      table_mcs_loss(ms_t, w, loss = loss_func),      table_mcs_loss(ms_t_out_m, w, loss = loss_func),      table_mcs_loss(ms_t_out_w, w, loss = loss_func))
  
  tab_loss <- round(rbind(loss_dgp_garch, loss_dgp_figarch, loss_dgp_dcs, loss_dgp_sv, loss_dgp_ms), 4)
  tab_mcs_loss <- round(rbind(loss_mcs_dgp_garch, loss_mcs_dgp_figarch, loss_mcs_dgp_dcs, loss_mcs_dgp_sv, loss_mcs_dgp_ms), 4)
  
  row.names(tab_loss) <- rep(c("GARCH-N", "GARCH-T", "FIGARCH-N", "FIGARCH-T", "DCS-N", "DCS-T", "SV-N", "SV-T", "MS-N", "MS-T"), 5)
  
  
  tab_latex <- tab_loss |> kbl(format = "latex", booktabs = TRUE, digits = 3, format.args = list(scientific = FALSE))
  
  for (j in 1:ncol(tab_loss)) {
    bg <- ifelse(tab_mcs_loss[, j] == 1, "gray!45", "white")
    tab_latex <- tab_latex |> column_spec(j + 1, background = bg)
  }
  return(tab_latex)
}



probability_regime_given_time_n <- function(p, q, sigma, r, Pt) {
  numA <- (1 - q) * dnorm(r, 0, sigma[2]) * (1 - Pt)
  numB <- p * dnorm(r, 0, sigma[1]) * Pt
  deno <- dnorm(r, 0, sigma[1]) * Pt + dnorm(r, 0, sigma[2]) * (1 - Pt)
  l <- numA / deno + numB / deno
  return(l)
}

probability_regime_given_time_t <- function(p, q, sigma, r, Pt, nu) {
  numA <- (1 - q) * sqrt(nu / (nu - 2)) / sigma[2] * dt(r * sqrt(nu / (nu - 2)) / sigma[2], nu) * (1 - Pt)
  numB <- p * sqrt(nu / (nu - 2)) / sigma[1] * dt(r * sqrt(nu / (nu - 2)) / sigma[1], nu) * Pt
  deno <- sqrt(nu / (nu - 2)) / sigma[1] * dt(r * sqrt(nu / (nu - 2)) / sigma[1], nu) * Pt +
    sqrt(nu / (nu - 2)) / sigma[2] * dt(r * sqrt(nu / (nu - 2)) / sigma[2], nu) * (1 - Pt)
  l <- numA / deno + numB / deno
  return(l)
}

msgarch_fit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  opt_methods <- c("BFGS", "Nelder-Mead", "CG", "SANN", "solnp")
  expr <- NULL
  while (is_error == TRUE && k < 6) {
    k <- k + 1
    expr <- tryCatch({
      if (k < 5) {
        fit_model <- FitML(spec, data, ctr = list(do.se = FALSE, do.plm = FALSE, OptimFUN = function(vPw, f_nll, spec, data, do.plm){
          out <- stats::optim(vPw, f_nll, spec = spec, data = data, do.plm = do.plm, method = opt_methods[k])}))
      } else {
        fit_model <- FitML(spec, data, ctr = list(do.se = FALSE, do.plm = TRUE, OptimFUN = function(vPw, f_nll, spec, data, do.plm) {
          fn_obj <- function(pars) {
            val <- f_nll(pars, spec = spec, data = data, do.plm = do.plm)
            return(sum(val))
          }
          out <- Rsolnp::solnp(pars = vPw, fun = fn_obj)
          out$value <- out$values[length(out$values)]
          out$convergence <- ifelse(out$convergence == 0, 0, 1)
          return(out)
        }))
      }
      TRUE
    }, error = function(e) {
      FALSE
    }, warning = function(cond) {
      FALSE
    })
    if (isTRUE(expr)) is_error <- FALSE
  }
  return(fit_model)
}


estimate_parameters_t <- function (data) {
  obj <- get_nll(data, model = "t", silent = TRUE, hessian = TRUE)
  fit <- stats::nlminb(obj$par, obj$fn, obj$gr,  lower = c(NULL, -5, NULL, NULL), upper = c(NULL, NULL, NULL,  5.52))
  rep <- suppressWarnings(TMB::sdreport(obj))
  k <- 0
  while (fit$convergence != 0 || any(is.nan(rep$sd)) || !rep$pdHess && k < 1000) {
    k <- k + 1
    obj <- get_nll(data, model = "t", silent = TRUE, hessian = TRUE)
    obj$par <- fit$par + runif(length(obj$par), -0.5, 0.5)
    fit <- stats::nlminb(obj$par, obj$fn, obj$gr, lower = c(NULL, -5 + runif(1, 0, 1) , NULL, NULL), upper = c(NULL, NULL, NULL, 5.52 + runif(1, -1, 1)))
    rep <- suppressWarnings(TMB::sdreport(obj))
  }
  
  if (fit$convergence != 0 || any(is.nan(rep$sd)) || !rep$pdHess) {
    stop("TMB optimization failed in estimate_parameters_n()")
  }
  
  opt <- list()
  class(opt) <- "stochvolTMB"
  opt$rep <- rep
  opt$obj <- obj
  opt$fit <- fit
  opt$nobs <- length(data)
  opt$model <- "t"
  opt$data <- data
  return(opt)
}

estimate_parameters_n <- function (data) {
  obj <- get_nll(data, model = "gaussian", silent = TRUE, hessian = TRUE)
  fit <- stats::nlminb(obj$par, obj$fn, obj$gr,  lower = c(NULL, -5, NULL))
  rep <- suppressWarnings(TMB::sdreport(obj))
  k <- 0
  while (fit$convergence != 0 || any(is.nan(rep$sd)) || !rep$pdHess && k < 1000) {
    k <- k + 1
    obj <- get_nll(data, model = "gaussian", silent = TRUE, hessian = TRUE)
    obj$par <- fit$par + runif(length(obj$par), -0.1, 0.1)
    fit <- stats::nlminb(obj$par, obj$fn, obj$gr, lower = c(NULL, -5, NULL), upper = c(NULL, NULL, NULL))
    rep <- suppressWarnings(TMB::sdreport(obj))
  }
  
  if (fit$convergence != 0 || any(is.nan(rep$sd)) || !rep$pdHess) {
    stop("TMB optimization failed in estimate_parameters_n()")
  }
  
  opt <- list()
  class(opt) <- "stochvolTMB"
  opt$rep <- rep
  opt$obj <- obj
  opt$fit <- fit
  opt$nobs <- length(data)
  opt$model <- "gaussian"
  opt$data <- data
  return(opt)
}

dcsn_fit <- function(data) {
  par_ini <- grid_dcsn(data)
  opt_methods <- c("BFGS", "Nelder-Mead", "CG", "SANN")
  is_error <- TRUE
  k <- 0
  params <- NULL
  while (is_error == TRUE && k < 4) {
    k <- k + 1
    expr <- tryCatch({
        fit <- optim(par = par_ini, dcsn_like, r = data, method = opt_methods[k])
        if (fit$convergence != 0) stop("No convergence")
        params <- fit$par
        TRUE
      }, error = function(e) {
        FALSE
      })
    if (isTRUE(expr)) is_error <- FALSE
  }
  return(params)
}

dcsn_lev_fit <- function(data) {
  par_ini <- grid_dcsn_lev(data)
  opt_methods <- c("BFGS", "Nelder-Mead", "CG", "SANN")
  is_error <- TRUE
  k <- 0
  params <- NULL
  while (is_error == TRUE && k < 4) {
    k <- k + 1
    expr <- tryCatch({
      fit <- optim(par = par_ini, dcsn_lev_like, r = data, method = opt_methods[k])
      if (fit$convergence != 0) stop("No convergence")
      params <- fit$par
      TRUE
    }, error = function(e) {
      FALSE
    })
    if (isTRUE(expr)) is_error <- FALSE
  }
  return(params)
}

vol_dcsn <- function(r, params) {
  n <- length(r) + 1
  lambda <- rep(0, n)
  lambda_star <- rep(0, n)
  lambda_star[1] <- 0
  lambda[1] <- params[1]
  for (i in 2:n) {
    u <- r[i - 1]^2 / exp(2 * lambda[i - 1]) - 1
    lambda_star[i] <- params[2] * lambda_star[i - 1] + params[3] * u
    #lambda_star[i] <- params[2] * lambda_star[i - 1] + params[3] * (r[i - 1]^2 / exp(2 * lambda[i - 1]) - 1)
    lambda[i] <- params[1] + lambda_star[i]
  }
  return(exp(lambda))
}

vol_dcsn_lev <- function(r, params) {
  n <- length(r) + 1
  lambda <- rep(0, n)
  lambda_star <- rep(0, n)
  lambda_star[1] <- 0
  lambda[1] <- params[1]
  for (i in 2:n) {
    u <- r[i - 1]^2 / exp(2 * lambda[i - 1]) - 1
    lambda_star[i] <- params[2] * lambda_star[i - 1] + params[3] * u + params[4] * sign(-r[i - 1]) * (u + 1)
    lambda[i] <- params[1] + lambda_star[i]
  }
  return(exp(lambda))
}

figarch_fit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  n <- length(data)
  opt_methods <- c("nloptr", "solnp", "gosolnp", "lbfgs", "nlminb", "hybrid")
  while (is_error && k < length(opt_methods)) {
    k <- k + 1
    is_error <- FALSE
    expr <- tryCatch({
      fit_model <- ugarchfit(spec, data, solver = opt_methods[k], fit.control = list(trunclag = min(n - 100, 1000)))
      warns <- character(0)
      fc <- withCallingHandlers(
        ugarchforecast(fit_model, n.ahead = 1),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")}
      )
        if (fit_model@fit$convergence != 0) stop("no convergence")
        if (any(!is.finite(fc@forecast$sigmaFor))) stop("NaN or Inf in forecast")
        if (any(grepl("Positivity Contraints", warns))) stop("Positivity Constraints NOT")
        TRUE
    }, error = function(e) {
        FALSE
      })
    if (!isTRUE(expr)) is_error <- TRUE
  }
  
  if (is_error) stop("All optimization methods failed to converge.")
  return(fit_model)
}


figarch_larger <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  n <- length(data)
  opt_methods <- c("nloptr", "solnp", "gosolnp", "lbfgs", "nlminb", "hybrid")
  while (is_error && k < length(opt_methods)) {
    k <- k + 1
    is_error <- FALSE
    expr <- tryCatch({
      fit_model <- ugarchfit(spec, data, solver = opt_methods[k], fit.control = list(trunclag = 2000))
      warns <- character(0)
      fc <- withCallingHandlers(
        ugarchforecast(fit_model, n.ahead = 1),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")}
      )
      if (fit_model@fit$convergence != 0) stop("no convergence")
      if (any(!is.finite(fc@forecast$sigmaFor))) stop("NaN or Inf in forecast")
      if (any(grepl("Positivity Contraints", warns))) stop("Positivity Constraints NOT")
      TRUE
    }, error = function(e) {
      FALSE
    })
    if (!isTRUE(expr)) is_error <- TRUE
  }
  
  if (is_error) stop("All optimization methods failed to converge.")
  return(fit_model)
}

garch_fit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  n <- length(data)
  opt_methods <- c("hybrid", "nloptr", "solnp", "gosolnp", "lbfgs", "nlminb")
  while (is_error && k < length(opt_methods)) {
    k <- k + 1
    is_error <- FALSE
    expr <- tryCatch({
      fit_model <- ugarchfit(spec, data, solver = opt_methods[k])
      warns <- character(0)
      fc <- withCallingHandlers(
        ugarchforecast(fit_model, n.ahead = 1),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")}
      )
      if (fit_model@fit$convergence != 0) stop("no convergence")
      if (any(!is.finite(fc@forecast$sigmaFor))) stop("NaN or Inf in forecast")
      if (any(grepl("Positivity Contraints", warns))) stop("Positivity Constraints NOT")
      TRUE
    }, error = function(e) {
      FALSE
    })
    if (!isTRUE(expr)) is_error <- TRUE
  }
  
  if (is_error) stop("All optimization methods failed to converge.")
  return(fit_model)
}




gas_fit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  opt_methods <- c("BFGS", "Nelder-Mead", "CG", "SANN")
  fit_model <- NULL
  while (is_error == TRUE && k < length(opt_methods)) {
    k <- k + 1
    tryCatch(
      expr <- {
        fit_model <- UniGASFit(spec, data, Compute.SE = FALSE, fn.optimizer = function(par0, data, GASSpec, FUN) {
          optimizer = optim(par0, FUN, data = data, GASSpec = GASSpec, method = opt_methods[k],
            control = list(trace = 0), hessian = FALSE)
          out = list(pars = optimizer$par,
            value = optimizer$value,
            hessian = optimizer$hessian,
            convergence = optimizer$convergence)
          return(out)
        })
        if (fit_model@Estimates$optimiser$convergence != 0) stop("No convergence")
        TRUE
      },
      error = function(e) {
        FALSE
      })
    if (isTRUE(expr)) is_error <- FALSE
  }
  return(fit_model)
}
  


library(quantreg)
VaR_VQR = function(r,VaR, alpha){
}


calibration_tests <- function(r_oos, var, es, s, alpha) {
  VaRBack <- BacktestVaR(r_oos, var, alpha = alpha, Lags = 4)
  aux_tests <- c(VaRBack$LRuc[2], VaRBack$LRcc[2], 
    cc_backtest(r_oos,  var, es, s, alpha  = alpha)$pvalue_twosided_general, 
    esr_backtest(r_oos, var, es, alpha  = alpha, B = 0, version = 1)$pvalue_twosided_asymptotic)
  return(sum(aux_tests > 0.05))
}

calibration_tests_describe <- function(r_oos, var, es, s, alpha) {
  VaRBack <- BacktestVaR(r_oos, var, alpha = alpha, Lags = 4)
  aux_tests <- round(c(as.numeric(VaRBack$LRuc[2]), as.numeric(VaRBack$LRcc[2]), VaRBack$DQ$pvalue,
    cc_backtest(r_oos,  var, es, s, alpha  = alpha)$pvalue_twosided_general, 
    esr_backtest(r_oos, var, es, alpha  = alpha, B = 0, version = 1)$pvalue_twosided_asymptotic), 3)
  aux <- c(ifelse(aux_tests[1] < 0.05, paste(aux_tests[1], "UC"), "NA"),
      ifelse(aux_tests[2] < 0.05, paste(aux_tests[2], "CC"), "NA"),
      ifelse(aux_tests[3] < 0.05, paste(aux_tests[3], "DQ"), "NA"),
      ifelse(aux_tests[4] < 0.05, paste(aux_tests[4], "CoC"), "NA"),
      ifelse(aux_tests[5] < 0.05, paste(aux_tests[5], "ESR"), "NA"))
  return(aux)
}


mcs_scoring_functions <- function(r_oos, var, es, alpha) {
  
  loss_ql  <- apply(var, 2, function(v) QL(matrix(v, ncol = 1), r_oos, alpha = alpha))
  loss_fzg <- sapply(seq_len(ncol(var)), function(i) FZG(matrix(var[, i], ncol = 1), matrix(es[, i],  ncol = 1), r_oos, alpha = alpha))
  loss_nz  <- sapply(seq_len(ncol(var)), function(i) NZ(matrix(var[, i], ncol = 1), matrix(es[, i],  ncol = 1), r_oos, alpha = alpha))
  loss_al  <- sapply(seq_len(ncol(var)), function(i) AL(matrix(var[, i], ncol = 1), matrix(es[, i],  ncol = 1), r_oos, alpha = alpha))
  
  table_mcs <- matrix(0, ncol = ncol(var), nrow = 4)
  table_mcs[1, estMCS.quick(loss_ql , test = "t.range", B = 10000, l = 21, alpha = 0.25)] <- 1
  table_mcs[2, estMCS.quick(loss_fzg , test = "t.range", B = 10000, l = 21, alpha = 0.25)] <- 1
  table_mcs[3, estMCS.quick(loss_nz , test = "t.range", B = 10000, l = 21, alpha = 0.25)] <- 1
  table_mcs[4, estMCS.quick(loss_al , test = "t.range", B = 10000, l = 21, alpha = 0.25)] <- 1
  
  return(apply(table_mcs, 2, sum))
}

