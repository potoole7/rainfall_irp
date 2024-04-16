#### Model fitting functions ####

# look at correlations between variables
plot_corr <- function(dat, title) {
  cors <- dat %>% 
    mutate(date = as.numeric(date)) %>% 
    select(-c(
      contains("min"), 
      contains("max"), 
      any_of(c("more_rain", "lat", "lon", "alt"))
    )) %>% 
    cor()
  
  cors %>% 
   ggcorrplot::ggcorrplot(
     outline.color = "white",
     ggtheme = theme_bw(),
     tl.cex = 8, 
     title = title,
     legend.title = "Correlation"
   ) %>% 
   print()
  return(cors)
}

# prefer fit_nested over fit_complex? https://api.rpubs.com/tomanderson_34/lrt
# Function to perform Likelihood ratio test (on negative log likelihood!)
lrtest <- function(fit_complex, fit_nested) {
  if (!"nllh" %in% names(fit_nested)) fit_nested <- fit_nested$fit
  teststat <- -2 * (-fit_nested$nllh + fit_complex$nllh)
  return(pchisq(teststat, df = 1, lower.tail = FALSE))
}

# also calculate AIC for a model (https://www.statology.org/aic-in-r/)
aic <- function(fit, k = NULL) {
  if (is.null(k)) k <- length(fit$mle) - 2
  # AIC = 2K – 2ln(L), K defaults to 2
  2 * (k + 2) - 2 * (- fit$nllh)
}

# Plot MLEs
# TODO: Improve!
plot_mles <- function(fit, predictors) {
  
  mles <- as.vector(fit$mle)
  if (is.null(predictors)) {
    names(mles) <- c("scale", "shape")
  # if stationary fit, mles are just for shape, scale varies over time
  } else {
    names(mles) <- c("intercept", predictors, "shape")
  }
  
  ses <- fit$se
  names(ses) <- names(mles)
  
  # plot 
  tibble(
    "name" = factor(names(mles), levels = names(mles)), 
    "mle"  = mles, 
    "se"   = ses
  ) %>% 
    ggplot(aes(x = name, y = mle)) + 
    geom_point() + 
    # should I be using se or se / 2 for error bounds around estimates?
    geom_errorbar(aes(ymin = mle - se, ymax = mle + se)) + 
    geom_hline(yintercept = 0, colour = "red") + 
    theme_bw()
}

# test adding one predictor and refitting
fit_models <- function(predictors, predict_mat_current, predict_mat) {
  lapply(seq_along(predictors), function(i) {
    ydat <- predict_mat[, i, drop = FALSE]
    sigl <- 1
    if (!is.null(predict_mat_current)) {
      ydat <- cbind(predict_mat_current, ydat)
      sigl <- ncol(predict_mat_current) + 1
    }
    ret <- ismev::gpd.fit(
      xdat      = dat$abs_diff,
      threshold = threshold,
      # ydat      = predict_mat[, i],
      ydat      = ydat,
      # sigl      = 1,
      sigl      = sigl,
      siglink   = exp,
      show      = FALSE
    )
    ret$pred <- predictors[-i]
    return(ret)
  })
}

# peform log likelihood test for each nested model of orig_fit
# Also compare using AIC
compare_models <- function(orig_fit, fit_lst, type = "lrt") {
  if (type == "lrt") {
    vapply(fit_lst, lrtest, orig_fit, FUN.VALUE = numeric(1))
  } else if (type == "aic") {
    vapply(fit_lst, aic, FUN.VALUE = numeric(1))
  }
}

# perform forward step (i.e. compare current model to those with 
# 1 additional parameter)
forward_step <- function(
    orig_fit, 
    all_predictors, 
    best_predictors = c(), 
    predict_mat, 
    best_pred_mat = NULL, 
    type = "lrt",
    sig = 0.05
) {
  
  # initialise best fit and associated predictors and matrix of predictors
  best_fit <- orig_fit
  
  # fit models with each predictor added
  fit_lst <- fit_models(all_predictors, best_pred_mat, predict_mat)
  # perform negative log likelihood test
  if ("fit" %in% names(best_fit)) {
    test_stats <- compare_models(best_fit$fit, fit_lst, type = type)
  } else test_stats <- compare_models(best_fit, fit_lst, type = type)
  # select best model
  if (type == "lrt") {
    best_mod <- which(test_stats < sig & test_stats == min(test_stats))
    # AIC
  } else {
    if ("fit" %in% names(best_fit)) {
      best_mod <- which.min(c(aic(best_fit$fit), test_stats))
    } else {
      best_mod <- which.min(c(aic(best_fit), test_stats))
    }
    # keep original model if identified as best model
    if (best_mod == 1) best_mod <- NULL else best_mod <- best_mod - 1
  }
  if (length(best_mod) > 0) {
    print(type)
    print(paste0(
      ifelse(type == "lrt", "Likelihood ratio test ", "AIC "),
      "selects model with predictor for `", 
      all_predictors[best_mod], 
      "`"
    ))
    best_fit <- fit_lst[[best_mod]]
    best_predictors <- c(best_predictors, predictors[best_mod])
    best_pred_mat <- predict_mat[, best_predictors]
  } else {
    # print("Likelihood ratio test selects original model")
    print(paste0(
      ifelse(type == "lrt", "Likelihood ratio test ", "AIC "),
      "selects model with original model"
    ))
    # TODO: Return table of likelihood ratio tests for each predictor
  }
  
  return(list(
    "fit"        = best_fit, 
    "predictors" = best_predictors, 
    "pred_mat"   = best_pred_mat,
    "test_stats" = test_stats
  ))
}

# find best models, using forward stepwise regression
# TODO: Return p-values for each step here
find_best_mod <- function(
  orig_fit, 
  all_predictors, 
  predict_mat, 
  type = "lrt", 
  sig = 0.05
) {
  best_fit <- orig_fit
  best_predictors <- c()
  best_pred_mat <- NULL
  fitted <- FALSE
  while(fitted == FALSE) {
    loop_fit <- forward_step(
      best_fit, 
      predictors, 
      best_predictors, 
      predict_mat, 
      best_pred_mat, 
      type = type, 
      sig = sig
    )
    # stopping criteria
    if (type == "lrt") {
      cond <- all(loop_fit$test_stats > sig)
    } else {
      if ("fit" %in% names(best_fit)) {
        aic <- aic(best_fit$fit)
      } else {
        aic <- aic(best_fit)
      }
      # 
      cond <- which.min(c(aic, loop_fit$test_stats)) == 1
    }
    if (cond) {
      fitted <- TRUE
    } else {
      best_fit <- loop_fit
      best_predictors <- best_fit$predictors
      best_pred_mat <- best_fit$pred_mat
    }
  }
  return(list(
    "best_fit"        = best_fit,
    "best_predictors" = best_predictors
  ))
}


#### Stolen functions ####

#' Threshold selection method for univariate extremes
#' https://github.com/conor-murphy4/automated_threshold_selection/tree/main

#' Threshold selection method for univariate extremes
#'
#' 'thresh_qq_metric' selects a constant threshold above which the data can be most closely modelled by a Generalised Pareto distribution.
#'
#' @author Conor Murphy
#'
#' @param data A numeric vector.
#' @param thresh A numeric vector of proposed thresholds to test.
#' @param k  A positive integer denoting the number of bootstraps.
#' @param m A positive integer denoting the number of equally-spaced probabilities at which to evaluate quantiles.
#'
#' @returns A list containing the chosen threshold, the parameters of the fitted GPD, the number of observations above the chosen thresholds and the metric values 'd' corresponding to each proposed threshold.
#'
#' @examples
#' set.seed(12345)
#' data_test1 <- rgpd(1000, shape = 0.1, scale = 0.5, mu = 1)
#' thresholds1 <- quantile(data_test1, seq(0, 0.95, by = 0.05))
#' (example1 <- thresh_qq_metric(data_test1, thresh = thresholds1))
#'
#' set.seed(11111)
#' test2 <- rgpd(10000, shape = 0.1, scale = 0.5)
#' u <- 1
#' cens_thr <- u * rbeta(length(test2), 1, 0.5)
#' keep <- test2 > cens_thr
#' data_test2 <- test2[keep]
#' thresholds2 <- quantile(data_test2, seq(0, 0.95, by = 0.05))
#' (example2 <- thresh_qq_metric(data_test2, thresh = thresholds2))
thresh_qq_metric <- function(data, thresh, k = 100, m = 500) {
  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.numeric(thresh)) stop("u to be tested needs to be a vector")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")

  meandistances <- xis <- sigmas <- num_excess <- numeric(length(thresh))
  for (i in 1:length(thresh)) {
    u <- thresh[i]
    excess <- data[data > u] - u
    num_excess[i] <- length(excess)
    if (num_excess[i] > 10) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL, z = excess, par = c(mle0, 0.1), control = list(fnscale = -1))
      xis[i] <- init.fit$par[[2]]
      sigmas[i] <- init.fit$par[[1]]
      distances <- numeric(k)
      for (j in 1:k) {
        X <- sample(excess, num_excess[i], replace = TRUE)
        mle <- mean(X)
        ifelse(xis[i] < 0, pars_init <- c(mle, 0.1), pars_init <- c(sigmas[i], xis[i]))
        gpd.fit <- optim(GPD_LL, z = X, par = pars_init, control = list(fnscale = -1))
        quants <- qgpd((1:m) / (m + 1), scale = gpd.fit$par[[1]], shape = gpd.fit$par[[2]])
        distances[j] <- (1 / m) * sum(abs(quantile(X, probs = (1:m) / (m + 1)) - quants))
      }
      meandistances[i] <- mean(distances)
    } else {
      meandistances[i] <- NA
    }
  }
  chosen_index <- which.min(meandistances)
  chosen_threshold <- thresh[chosen_index]
  xi <- xis[chosen_index]
  sigma <- sigmas[chosen_index]
  len <- num_excess[chosen_index]
  result <- list(thresh = chosen_threshold, par = c(sigma, xi), num_excess = len, dists = meandistances)
  return(result)
}

# =====================================================================
# Functions for Generalised Pareto Distribution.
# Added option to use nu parameterisation
# checks that param values are valid
# =====================================================================
# pgpd
# qgpd
# dgpd
# rgpd

# GPD_LL
# =====================================================================
#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p<mu.
#'
#' @author Zak Varty
#'
#' @param q vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @param skip_checks logical. Speed up evaluation by skipping checks on inputs? (Beware!)
#' @return Probability of the GPD X<=q
#' @importFrom stats pexp
#' @examples
#' pgpd(q = c(-1, 1.5, 3), shape = 1, scale = 1)
#' pgpd(q = 1.5, shape = c(0, -1), scale = c(0.1, 1))
#' @export

pgpd <- function(q, shape, scale = NULL, nu = NULL, mu = 0, skip_checks = FALSE) {
  if (!skip_checks) {
    # one and only one of {nu, scale} may be specified
    if (is.null(scale) & is.null(nu)) {
      stop("Define one of the parameters nu or scale.")
    }
    if (!is.null(scale) & !is.null(nu)) {
      stop("Define only one of the parameters nu and scale.")
    }
    # Calculate scale from nu if required
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
      if (any(scale <= 0)) {
        stop("Implied scale parameter(s) must be positive.")
      }
    }
    # Check that scale value(s) are positive
    if (any(scale <= 0)) {
      stop("Scale parameter(s) must be positive.")
    }

    # Ensure q, scale, shape and mu are of same length.
    if (length(scale) == 1 & length(q) > 1) {
      scale <- rep(scale, length(q))
    }
    if (length(shape) == 1 & length(q) > 1) {
      shape <- rep(shape, length(q))
    }
    if (length(mu) == 1 & length(q) > 1) {
      mu <- rep(mu, length(q))
    }
  } else {
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
    }
  }
  # calculate probabilities
  p <- (1 - (1 + (shape * (q - mu)) / scale)^(-1 / shape))
  # correct probabilities below mu or above upper end point
  p[q < mu] <- 0
  p[(shape < 0) & (q >= (mu - scale / shape))] <- 1

  # correct probabilities where xi = 0
  if (any(abs(shape) < 1e-10)) {
    # ex <- which(shape ==0)
    ex <- which(abs(shape) < 1e-10)
    p[ex] <- pexp(q = q[ex] - mu[ex], rate = 1 / scale[ex])
  }

  return(p)
}

#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p is not a valid
#' probability.
#'
#' @author Zak Varty
#'
#' @param p vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @return Probability of the GPD X<=x
#' @examples
#' qgpd(p = 0.5, shape = 0.5, scale = 0.5)
#' \dontrun{
#' qgpd(p = -0.1, shape = 0, scale = 1, mu = 0.1)
#' }
#' @export
qgpd <- function(p, shape, scale = NULL, nu = NULL, mu = 0) {
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop("Define one of the parameters nu or scale.")
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop("Define only one of the parameters nu and scale.")
  }

  # Probabilities must all be positive
  if (!all((p >= 0) & (p <= 1))) {
    stop("Probabilities p must be in the range [0,1].")
  }

  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop("Implied scale parameter(s) must be positive.")
    }
  }

  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop("Scale parameter(s) must be positive.")
  }
  # Ensure p, scale, shape and mu are of same length.
  if (length(scale) == 1 & length(p) > 1) {
    scale <- rep(scale, length(p))
  }
  if (length(shape) == 1 & length(p) > 1) {
    shape <- rep(shape, length(p))
  }
  if (length(mu) == 1 & length(p) > 1) {
    mu <- rep(mu, length(p))
  }

  # calculate quantiles
  q <- mu + (scale / shape) * ((1 - p)^(-shape) - 1)

  # correct quantiles where xi = 0
  # ex <- which(shape ==0)
  if (any(abs(shape) < 1e-10)) {
    ex <- which(abs(shape) < 1e-10)
    q[ex] <- mu[ex] + stats::qexp(p = p[ex], rate = 1 / scale[ex])
  }
  return(q)
}



#' Generalised Pareto Distribution
#'
#' Density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or x is
#' outside of the domain of the given distribution.
#'
#' @author Zak Varty
#'
#' @param x vector of values as which to evaluate density.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter
#' @param mu  location parameter
#' @param log  locical. Return log
#' @return density of the GPD at x
#' @examples
#' dgpd(x = c(-1, 0.5, 1, 1.9, 5), shape = -0.5, scale = 1)
#' @export
#'
dgpd <- function(x, shape, scale = NULL, nu = NULL, mu = 0, log = FALSE) {
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop("Define one of the parameters nu or scale.")
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop("Define only one of the parameters nu and scale.")
  }

  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop("Implied scale parameter(s) must be positive.")
    }
  }

  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop("Scale parameter(s) must be positive.")
  }
  # Ensure x, scale, shape and mu are of same length.
  if (length(scale) == 1 & length(x) > 1) {
    scale <- rep(scale, length(x))
  }
  if (length(shape) == 1 & length(x) > 1) {
    shape <- rep(shape, length(x))
  }
  if (length(mu) == 1 & length(x) > 1) {
    mu <- rep(mu, length(x))
  }

  if (log == FALSE) {
    out <- (scale^(-1)) * pmax((1 + shape * (x - mu) / scale), 0)^((-1 / shape) - 1)
    # amend values below threshold
    out[which(x < mu)] <- 0
    # amend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale / shape)))] <- 0
    # amend values where xi = 0 (if they exist)
    if (any(abs(shape < 1e-10))) {
      ex <- which(abs(shape) < 1e-10)
      out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1 / scale[ex])
    }
  } else {
    out <- -log(scale) + ((-1 / shape) - 1) * log(pmax((1 + shape * (x - mu) / scale), 0))
    # amend values below threshold
    out[which(x < mu)] <- -Inf
    # amend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale / shape)))] <- -Inf
    # amend values where xi = 0 (if they exist)
    if (any(abs(shape) < 1e-10)) {
      ex <- which(abs(shape) < 1e-10)
      out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1 / scale[ex], log = TRUE)
    }
  }
  return(out)
}

#' Generalised Pareto Distribution
#'
#' Sample the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0.
#'
#' @author Zak Varty
#'
#' @param n sample size.
#' @param shape shape parameter (xi).
#' @param scale scale parameter (sigma).
#' @param nu  alternative scale parameter.
#' @param mu  location parameter.
#' @return Random sample from generalised pareto distirbution.
#'
#' @examples
#' rgpd(n = 100, shape = 0, scale = 1:100)
#' @export
rgpd <- function(n, shape, scale = NULL, nu = NULL, mu = 0) {
  ## Input checks
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop("Define one of the parameters nu or scale.")
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop("Define only one of the parameters nu and scale.")
  }
  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop("Implied scale parameter(s) must be positive.")
    }
  }
  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop("Scale parameter(s) must be positive.")
  }
  # Ensure q, scale, shape and mu are of same length.
  if ((length(scale) == 1) & (n > 1)) {
    scale <- rep(scale, n)
  }
  if ((length(shape) == 1) & (n > 1)) {
    shape <- rep(shape, n)
  }
  if ((length(mu) == 1) & (n > 1)) {
    mu <- rep(mu, n)
  }

  # simulate sample
  sample <- mu + (scale / shape) * ((1 - stats::runif(n))^(-shape) - 1)
  # correct sample values where xi = 0
  # ex <- which(shape ==0)
  if (any(abs(shape) < 1e-10)) {
    ex <- which(abs(shape) < 1e-10)
    sample[ex] <- mu[ex] +
      stats::rexp(n = length(ex), rate = 1 / scale[ex])
  }
  return(sample)
}

# rgpd_rd <- function(n, sig, xi, mu, to_nearest, mu_latent = NULL){
#  if(is.null(mu_latent)) mu_latent = mu - 0.5 * to_nearest
#  x <- rgpd(n = n, shape = xi, scale = sig, mu = mu_latent)
#  x <- round_to_nearest(x, to_nearest)
#  return(x)
# }

#' evalue probability mass function of rounded generalised Pareto distribution
#'
#' @author Zak Varty
#'
#' @param x Vector values at which to evaluate mass function
#' @param u Vector of latent threshold values
#' @param sig_u Vector of latent scale parameters (for exceedances of u)
#' @param xi Latent shape parameter
#' @param to_nearest Level of rounding
#'
#' @return pmf evaluated at x. NOTE: does not check validity of x values.
#'
#' @examples
#' dgpd_rd(x = seq(0.1, 1.5, by = 0.1), to_nearest = 0.1, u = 0, sig_u = 1, xi = 0)
#' dgpd_rd(x = seq(0.1, 1.5, by = 0.1), to_nearest = 0.1, u = seq(0, 1.4, by = 0.1), sig_u = 1, xi = 0)
#' # CAUTION:
#' gpd_rd(x = 0.15, to_nearest = 0.1, u = 0, sig_u = 1, xi = 0)
dgpd_rd <- function(x, u, sig_u, xi, to_nearest) {
  # If (Y_i - u_i | Y_i > u_i) ~ GPD(sig_i, xi)
  # then Z_i  = [(Y_i - u_i)/sig_i | Y_i - u_i > 0] ~ GPD(1, xi)

  # range of z values that lead to observing x
  x_low <- pmax(x - to_nearest / 2, u)
  z_low <- (x_low - u) / sig_u
  z_high <- (x - u + to_nearest / 2) / sig_u

  # calculate probability of z in that range
  p_high <- pgpd(q = z_high, scale = 1, shape = xi, mu = 0)
  p_low <- pgpd(q = z_low, scale = 1, shape = xi, mu = 0)
  p <- p_high - p_low

  return(p)
}


#' Generalised Pareto log-likelihood
#'
#' @author Conor Murphy
#'
#' @param par A numeric vector of parameter values of length 2.
#' @param z A numeric vector of excesses of some threshold.
#'
#' @returns A numeric value of the log-likeihood.
#'
#' @examples
#' test1 <- rgpd(1000, shape = 0.1, scale = 0.5, mu = 1)
#' excess <- test1[test1 > 1.5] - 1.5
#' GPD_LL(par = c(1, 0.4), z = excess)
GPD_LL <- function(par, z) {
  sigma <- par[1]
  xi <- par[2]
  if (sigma > 0) {
    if (abs(xi) < 1e-10) {
      return(-length(z) * log(sigma) - ((1 / sigma) * sum(z)))
    } else {
      if (all(1 + (xi * z) / sigma > 0)) {
        return(-(length(z) * log(sigma)) - ((1 + 1 / xi) * (sum(log(1 + (xi * z) / sigma)))))
      } else {
        return(-1e6)
      }
    }
  } else {
    return(-1e7)
  }
}

transform_to_exp <- function(y, sig, xi) {
  std_exp <- (1 / xi) * log(1 + xi * (y / sig))
  return(std_exp)
}