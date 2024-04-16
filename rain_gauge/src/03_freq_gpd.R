#### Fit (Frequentist) GPD to absolute differences #####

# Problems: 
# May not be estimating

# TODO: Could also provide predictors for shape parameter!?
# TODO: Also identify best model via AIC (currently working on this)
# TODO: Investigate quadratic relationship of temperature with Sam, compare predictions 
# TODO: Possible seasonality to temporal effect, which could be added?
# TODO: Make predictions to assess performance! R^2 (not valid), train/test, etc 
# (look into assessing EVT models)
# TODO: Perform Bayesian analysis with `brms`, make covariate effects more 
# interpretable

# TODO: Learn how to interpret profile log-likelihoods

# TODO: How to model location parameter with GPD? 

# TODOs for next meeting
# TODO: Make predictions using fitted model! Perhaps use these to perform 
# cross validation
# TODO: Produce return level plot for dependent model, using MCMC?
# TODO: Somehow simulate additional data? Just sample from GPD? Using same 
# environmental variables?

# TODO: Preprocess sonde data as per Anakin's paper, to get weather
# Compare to weather from gauge data

#### libs ####

library(tidyr)
library(lubridate)
library(ggplot2)
library(aod) # For Wald test
# EVT packages & functions
library(evd)
library(ismev)
devtools::load_all("~/forks/mev")
# https://github.com/conor-murphy4/automated_threshold_selection
source("rain_gauge/src/functions.R")

# mev also has filter function, annoyingly!
library(dplyr)

#### Load Data ####

# TODO: Define data directory in metadata section
gauge <- readr::read_csv("rain_gauge/data/rain_gauge_daily.csv.gz")
sat <- readr::read_csv("rain_gauge/data/sat_dat_daily.csv.gz")

# sat_22_23 <- readr::read_csv(
#   "rain_gauge/data/2022_2023/sat_dat_daily_2022_2023.csv.gz"
# )
# sat <- sat %>% 
#   bind_rows(sat_22_23) %>% 
#   arrange(date)

# sonde <- readr::read_csv(
#   "rain_gauge/data/sonde_daily.csv.gz"
# )
# sonde_22_23 <- readr::read_csv(
#   "rain_gauge/data/2022_2023/sonde_daily_2022_2023.csv.gz"
# )
# sonde <- sonde %>% 
#   bind_rows(sonde_22_23) %>% 
#   arrange(date)

# merge data
dat <- full_join(
  x = gauge,
  y = sat %>% 
    select(date, sat_precip = precip)
) %>% 
  # calculate difference in rainfall
  mutate(
    diff = tbrg_precip_total - sat_precip,
    more_rain = case_when(
      diff < 0 ~ "satellite",
      diff > 0 ~ "gauge", 
      TRUE     ~ "neither"
    ), 
    abs_diff = abs(diff)
  ) %>% 
  # Temp: Remove dates which don't have data from both sources
  dplyr::filter(!is.na(abs_diff))

# TEMP: Remove all 0s
# dat <- dat %>% 
#   dplyr::filter(abs_diff != 0)

#### Parameter Stability Plots & Threshold Selection ####

# Mean residual life plot
par(mfrow = c(1, 1))
mrlplot(dat$abs_diff) # any threshold between > 0 and ~12 seems reasonable

# Threshold Choice Plot
par(mfrow = c(1, 2))
tcplot(dat$abs_diff, tlim = c(2, 8)) # max points with functioning uncertainty
tcplot(dat$abs_diff, tlim = c(0, 15), std.err = FALSE)
par(mfrow = c(1, 1))
# looks like anything over 6 is a good choice, so 7 seems reasonable

# thresholds <- quantile(dat$abs_diff, seq(0, 0.99, by = 0.01))
thresholds <- quantile(dat$abs_diff[dat$abs_diff != 0], seq(0, 0.99, by = 0.01))
# thcan <- quantile(dat$abs_diff, seq(0.8, 0.99, length.out = 25))
# thresh <- seq(0.8, 0.99, length.out = 25)

qq_thresholds <- thresh_qq_metric(dat$abs_diff, thresh = thresholds)

par(mfrow = c(1, 1))
plot(thresholds, qq_thresholds$dists, xlab="Threshold", ylab="Metric value")
abline(v = qq_thresholds$thresh, col="red", lwd = 2)
# gives choice of 1.754285 (or 16.85108 when 0s removed)
# number of exceedances: 47 (or 11 when 0s removed)

# thcan <- quantile(dat$abs_diff, seq(0, 0.99, length.out = 25))
# tstab.gpd(xdat = dat$abs_diff, thresh = thcan, method = "profile")
tstab.gpd(xdat = dat$abs_diff, thresh = thresholds, method = "profile")


# final choice of threshold
threshold <- qq_thresholds$thresh[[1]]
# TEMP: change to something higher
threshold <- 4 # better than above!!
# threshold <- 7
# threshold <- 6

# choose threshold based on repeated fitting
par(mfrow = c(1, 1))
# gpd.fitrange(dat$abs_diff, 1, 10, nint = 100)
gpd.fitrange(dat$abs_diff, 0, 6, nint = 100)
# again, seems to be stable for thresholds over 6
# However, we get NAs for shape when threshold is > 6!

# data above threshold
dat_remain <- dat %>% 
  dplyr::filter(abs_diff > threshold)

# plot absolute difference vs predictors
dat_remain %>% 
  select(
    abs_diff, c(date, temp_mean, rh_mean, atmos_pressure, wspd_vec_mean)
  ) %>% 
  mutate(date = as.numeric(date)) %>% 
  tidyr::pivot_longer(date:wspd_vec_mean) %>% 
  ggplot(aes(x = value, y = abs_diff)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~ name, scales = "free")
# seems to be some non-linear relationships going on here! 
# Might be worth exploring

# look at correlations between variables
plot_corr <- function(dat, title) {
  cors <- dat %>% 
    mutate(date = as.numeric(date)) %>% 
    select(-c(contains("min"), contains("max"), more_rain, lat, lon, alt)) %>% 
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

cors <- plot_corr(dat, "No Cutoff")
cors_remain <- plot_corr(dat_remain, paste0("Cutoff = ", round(threshold, 3)))

# possible predictors?
# predictors <- cors_remain[nrow(cors_remain), ]
# predictors <- names(predictors[abs(predictors) > 0.25])
# predictors <- predictors[
#   !grepl("precip|diff", predictors, fixed = FALSE)
# ]
# predictors # rh_mean, atmos_pressure, would like to have more!!

#### Fit simple generalised Pareto Distribution ####

fit1 <- ismev::gpd.fit(dat$abs_diff, threshold = threshold) # gives same nll

# Diagnostic plots
gpd.diag(fit1)
# Comments: QQ plot crap, definitely evidence of temporal dependence, 
# unsurprisingly
# return level plot is really weird, need to show to Christian!
# Also, estimate for scale parameter is very high (6.6258644), 
# perhaps we should also consider modelling it as temporally dependent

# profile log-likelihood for shape
gpd.profxi(fit1, xlow = floor(fit1$mle[1] - 2), xup = ceiling(fit1$mle[1] + 2))
# how to interpret? Definitely does not look like a "bell-shaped" distribution!

# profile log-likelihood for return level
gpd.prof(fit1, m = 10, xlow = 0, xup = max(fit1$data))
# how to interpret? Same comment!


#### GPD with non-stationary location #### 

# http://www.mas.ncl.ac.uk/~nlf8/teaching/mas8306/notes/chapter5.pdf

# prefer fit_nested over fit_complex? https://api.rpubs.com/tomanderson_34/lrt
# Function to perform Likelihood ratio test (on negative log likelihood!)
lrtest <- function(fit_complex, fit_nested) {
  teststat <- -2 * (-fit_nested$nllh + fit_complex$nllh)
  return(pchisq(teststat, df = 1, lower.tail = FALSE))
}

# also calculate AIC for a model (https://www.statology.org/aic-in-r/)
aic <- function(nllh, k) {
  # AIC = 2K – 2ln(L), K defaults to 2
  2 * (k + 2) - 2 * (- nllh)
}

# matrix of predictors, including date
# http://www.mas.ncl.ac.uk/~nlf8/teaching/mas8306/notes/chapter5.pdf
predict_mat <- dat %>% 
  mutate(
    date = as.numeric(date), 
    date = date + 1 - min(date) # start from 1
  ) %>% 
  select(date, temp_mean, rh_mean, atmos_pressure, wspd_vec_mean) %>% 
  as.matrix() %>% 
  # centre and scale, as per recommendation in `help(gpd.fit)`
  scale()

fit2 <- ismev::gpd.fit(
  dat$abs_diff, 
  threshold = threshold, 
  ydat = predict_mat,
  # model scale parameter with covariates
  sigl = c(1:ncol(predict_mat)),
  siglink = exp
)
# diagnostic plots
gpd.diag(fit2) # QQ plot seems to be much better!

# compare fit2 and fit1
# likelihood ratio test:
lrtest(fit2, fit1) 
# p-value << 0.05, reject null hypothesis to use nested model
# AIC
aic(fit2$nllh, k = ncol(predict_mat)) # 304.219
aic(fit1$nllh, k = 0) # 311.2635
# AIC lower for complete fit, prefer!

# Plot MLEs
# TODO: Improve plot!
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

# full set of predictors
predictors <- colnames(predict_mat)

# plot_mles(fit1, NULL)
# plot_mles(fit2, predictors)


#### Find best model ####

# Need to fit every nested model, to see how they compare
# TODO: Fix this, currently wrong!
best_fit <- fit2 # for now
best_fit$predictors <- predictors

# Test removing one predictor and refitting
fit_models <- function(predictors, predict_mat) {
  lapply(seq_along(predictors), function(i) {
    ret <- ismev::gpd.fit(
      xdat      = dat$abs_diff, 
      threshold = threshold, 
      ydat      = predict_mat[, -i],
      sigl      = c(1:(ncol(predict_mat) - 1)), 
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

# find and return best model, selected by likelihood ratio test
find_best_mod <- function(orig_fit, predictors, predict_mat, type = "lrt") {
  # initialise best fit and associated predictors and matrix of predictors
  best_fit <- orig_fit
  best_predictors <- predictors
  best_pred_mat <- predict_mat
  
  # fit models with each predictor removed
  for (i in seq_along(predictors)) {
    # fit models, removing each variable
    fit_lst <- fit_models(best_predictors, best_pred_mat)
    # perform negative log likelihood test
    lrtests <- compare_models(best_fit, fit_lst, type = type)
    # select best model
    best_mod <- which(lrtests < 0.05 & lrtests == min(lrtests))
    if (length(best_mod) > 0) {
      print(paste0(
        "Likelihood ratio test selects model without predictor for `", 
        best_predictors[best_mod], 
        "`"
      ))
      best_fit <- fit_lst[[best_mod]]
      best_predictors <- best_predictors[-best_mod]
      best_pred_mat <- best_pred_mat[, best_predictors]
    } else {
      # print("Likelihood ratio test selects original model")
      # TODO: Return table of likelihood ratio tests for each predictor
      break
    }
  }
  return(list(
    "fit" = best_fit, 
    "predictors"   = best_predictors, 
    "lrtest_pvals" = lrtests
  ))
}

# find best models
# best_fit <- find_best_mod(fit2, predictors, predict_mat)
# best_

# plot MLEs with standard errors
# plot_mles(best_fit$fit, best_fit$predictors)
plot_mles(best_fit, best_fit$predictors)

# gpd.diag(best_fit$fit) # why no return level plots?
gpd.diag(best_fit) # why no return level plots?

#### Predictions ####

# covariates for extreme observations
predict_mat_remain <- predict_mat[
  which(dat$abs_diff > threshold), 
]

# calculate fitted line for scale over time and with other covariates
# mles <- fit_best$mle
mles <- best_fit$mle
n <- ncol(predict_mat)

# TODO: can functionalise this recursively
sigma_t <- exp(
  mles[1] + 
    mles[2] * scale(predict_mat_remain[, 1]) + 
    mles[3] * scale(predict_mat_remain[, 2]) + 
    mles[4] * scale(predict_mat_remain[, 2]) + 
    mles[5] * scale(predict_mat_remain[, 3]) + 
    mles[6] * scale(predict_mat_remain[, 4]) + 
    mles[7] * scale(predict_mat_remain[, 5]) 
)

# make predictions
predictions <- rgpd(
  nrow(dat_remain), shape = last(best_fit$mle), scale = sigma_t
)

# Add predictions to data
dat_remain_pred <- dat_remain %>% 
  mutate(predictions = predictions)
  
# plot
dat_remain_pred %>% 
  select(date, abs_diff, predictions) %>% 
  pivot_longer(abs_diff:predictions) %>% 
  ggplot(aes(x = date, y = value, colour = name)) + 
  ggtitle(paste0("Threshold = ", threshold)) + 
  # geom_line() + 
  geom_point() + 
  theme_bw()


  
# R^2 only maxes sense for constant variance, 
# for GPD variance varies quite a bit!!
# Use WAIC, AIC, etc, likelihood ratio test
# Q-Q plot
# Assessing predictions: Use return levels, for multivariate we need to 
# specify covariates and sample from them using MCMC (look into this!)
# But I'm only looking at one location with non-stationary extremes, so this 
# should be fine http://www.mas.ncl.ac.uk/~nlf8/teaching/mas8391/background/chapter5.pdf
# Also sample from parameter values
