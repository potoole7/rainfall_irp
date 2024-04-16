#### Fit (Frequentist) GPD to absolute differences #####

# Fitting GPD with full dataset, rather than just 2022-2023
# Fitting with Sonde data, rather than weather data for 22-23 from gauge

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
# gauge <- readr::read_csv("rain_gauge/data/rain_gauge_daily.csv.gz")
gauge <- readr::read_csv("rain_gauge/data/rain_gauge_daily_full.csv.gz") %>% 
  rename(tbrg_precip_total = rain)
sat <- readr::read_csv("rain_gauge/data/sat_dat_daily.csv.gz")

# TEMP: Use 22-23 weather covariates for all years
# weather_covars <- readr::read_csv("rain_gauge/data/rain_covars.csv.gz")
# weather_covars <- readr::read_csv("rain_gauge/data/")
weather_covars <- readr::read_csv(
  "rain_gauge/data/sonde_daily_no_height.csv.gz"
) %>% 
  bind_rows(readr::read_csv(
    "rain_gauge/data/sonde_daily_no_height_22_23.csv.gz"
  )) %>% 
  select(-c(matches("height"), matches("bar_pres"))) %>% 
  # rename to match previous
  # rename(temp_mean = temp, rh_mean = rh, wspd_vec_mean = wspd) %>% 
  identity()

# quick plot of weather
# weather_covars %>% 
#   pivot_longer(rh:wspd) %>% 
#   ggplot(aes(x = date, y = value, colour = name)) + 
#   geom_line() + 
#   facet_wrap(~ name, scales = "free")
  
# gauge <- gauge %>% 
#   left_join(weather_covars)
  
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
  dplyr::filter(!is.na(abs_diff)) %>% 
  # join in weather covariates from sonde data
  left_join(weather_covars)

# couple days with NAs for weather! Remove
dat <- dplyr::filter(dat, !is.na(wspd))

# TEST: Investigate and remove two very large differences
dat %>% 
  dplyr::filter(abs_diff > 90)
# interestingly, gauge picks up rain not found by satellite, remove, note (?)
dat <- dat %>% 
  dplyr::filter(abs_diff < 90)

# save data for future analyses
readr::write_csv(dat, file = "rain_gauge/data/joined_rain_data.csv.gz")

#### Parameter Stability Plots & Threshold Selection ####

# plot abs differences themselves
dat %>% 
  dplyr::filter(abs_diff > 0) %>% 
  mutate(more_rain = ifelse(more_rain == "satellite", "radar", more_rain)) %>% 
  ggplot(aes(x = date, y = abs_diff, colour = more_rain)) + 
  geom_point()

# Mean residual life plot
par(mfrow = c(1, 1))
mrlplot(dat$abs_diff, tlim = c(0, 15)) 
# seem approx linear between ~3 and 7, quite unclear though

# Threshold Choice Plot
par(mfrow = c(1, 2))
tcplot(dat$abs_diff, tlim = c(0, 10)) # max points with functioning uncertainty
# uncertainty high, but seem to be roughly stable for threshold over 3 (?)
par(mfrow = c(1, 1))

# remove 0s, or keep?
thresholds <- quantile(dat$abs_diff, seq(0, 0.99, by = 0.01))
# thresholds <- quantile(dat$abs_diff[dat$abs_diff != 0], seq(0, 0.99, by = 0.01))
# thcan <- quantile(dat$abs_diff, seq(0.8, 0.99, length.out = 25))
# thresh <- seq(0.8, 0.99, length.out = 25)

# thcan <- quantile(dat$abs_diff, seq(0, 0.99, length.out = 25))
# tstab.gpd(xdat = dat$abs_diff, thresh = thcan, method = "profile")
par(mfrow = c(1, 2))
tstab.gpd(xdat = dat$abs_diff, thresh = thresholds, method = "profile")
par(mfrow = c(1, 1))
# prettier than above, choosing threshold of roughly 3 as well

set.seed(123)
qq_thresholds <- thresh_qq_metric(dat$abs_diff, thresh = thresholds)

par(mfrow = c(1, 1))
plot(thresholds, qq_thresholds$dists, xlab="Threshold", ylab="Metric value")
abline(v = qq_thresholds$thresh, col="red", lwd = 2)
print(paste("threshold:", round(qq_thresholds$thresh, 3)))
print(paste0(
  "number of exceedences: ", 
  nrow(dplyr::filter(dat, abs_diff > qq_thresholds$thresh)), 
  "/", 
  nrow(dat)
))

# final choice of threshold
threshold <- qq_thresholds$thresh[[1]]
# test different thresholds
# threshold <- 10
# threshold <- 5

# Alternative: choose threshold based on repeated fitting
par(mfrow = c(1, 1))
gpd.fitrange(dat$abs_diff, 0, 5, nint = 100)
# similar results to above, seems to become constant somewhere between 
# threshold of 2 and 3
 
# data above threshold
dat_remain <- dat %>% 
  dplyr::filter(abs_diff >= threshold)

# plot absolute difference vs predictors
dat_remain %>% 
  select(
    # abs_diff, c(date, temp, rh, atmos_pressure, wspd)
    abs_diff, c(date, temp, rh, wspd)
  ) %>% 
  mutate(date = as.numeric(date)) %>% 
  tidyr::pivot_longer(date:wspd) %>% 
  ggplot(aes(x = value, y = abs_diff)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~ name, scales = "free")

# Seems to be more exceedances for higher temperatures and certain 
# wind speeds
# rh doesn't seem to have much of an effect?
dat_plot <- dat %>% 
  dplyr::filter(!is.na(wdir)) %>% 
  select(-c(tbrg_precip_total, sat_precip))

cors <- plot_corr(dat_plot, "No Cutoff")
cors_remain <- plot_corr(select(dat_remain, -contains("precip")), paste0("Cutoff = ", round(threshold, 3)))

# possible predictors?
# predictors <- cors_remain[nrow(cors_remain), ]
# predictors <- names(predictors[abs(predictors) > 0.25])
# predictors <- predictors[
#   !grepl("precip|diff", predictors, fixed = FALSE)
# ]
# predictors # rh, atmos_pressure, would like to have more!!

#### Fit simple generalised Pareto Distribution ####

fit1 <- ismev::gpd.fit(dat$abs_diff, threshold = threshold) # gives same nll

# Diagnostic plots
gpd.diag(fit1)
# Prob plot good plot
# QQ plot generally good!
# return level plot seems to report very few occurences over the next 10 years

# profile log-likelihood for shape
gpd.profxi(fit1, xlow = 0, 3)
# how to interpret? 

# profile log-likelihood for return level
gpd.prof(fit1, m = 10, xlow = 0, xup = 20)
# how to interpret? Weird looking!


#### GPD with non-stationary location #### 

# http://www.mas.ncl.ac.uk/~nlf8/teaching/mas8306/notes/chapter5.pdf

# matrix of predictors, including date
# http://www.mas.ncl.ac.uk/~nlf8/teaching/mas8306/notes/chapter5.pdf
predict_mat <- dat %>% 
  mutate(
    date = as.numeric(date), 
    date = date + 1 - min(date) # start from 1
  ) %>% 
  # select(date, temp, rh, atmos_pressure, wspd) %>% 
  select(date, temp, rh, wspd) %>% 
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
gpd.diag(fit2) # QQ plot seems to be slightly better, off for two largest points
# How can I get return level plots for this?

# compare fit2 and fit1
# likelihood ratio test:
lrtest(fit2, fit1) 
# p-value < 0.001, reject null hypothesis to use nested model => use covariates!
# AIC
# complete model - nested model
aic(fit2, k = ncol(predict_mat)) - aic(fit1, k = 0) < 0
# AIC lower for complete fit, prefer!

# full set of predictors
predictors <- colnames(predict_mat)

# plot mles
plot_mles(fit1, NULL)
plot_mles(fit2, predictors) # positive, significant effect for all, especially temperature

#### Find best model ####

# Need to fit every nested model, to see how they compare
# TODO: Fix this, currently wrong! Do forward, rather than backwards stepwise

# initialise best model, using forward stepwise regression!
best_fit <- fit1

# find best model via lrt
best_fit_lrt <- find_best_mod(
    fit1, 
    predictors, 
    predict_mat, 
    sig = 0.05
)

# find best model according to AIC
# TODO: Fix for AIC
best_fit_aic <- find_best_mod(
    fit1, 
    predictors, 
    predict_mat, 
    type = "aic",
    sig = 0.05
)

# best fit is the same!
if (best_fit_lrt$best_fit$predictors == best_fit_lrt$best_fit$predictors) {
  print("Models selected via LRT and AIC are the same")
} 
best_fit <- best_fit_lrt

predictors_orig <- predictors
predictors <- best_fit$best_predictors
best_fit <- best_fit$best_fit$fit

# temp: 
best_fit <- fit2
predictors <- predictors_orig

# plot MLEs with standard errors
plot_mles(best_fit, predictors)
# plot_mles(best_fit, best_fit$predictors)

gpd.diag(best_fit) # why no return level plots?


#### Predictions ####

# covariates for extreme observations
predict_mat_remain <- predict_mat[
  which(dat$abs_diff >= threshold), 
]

# calculate fitted line for scale over time and with other covariates
# mles <- fit_best$mle
# mles <- best_fit$mle
mles <- best_fit$mle
n <- ncol(predict_mat)

# matrix of final predictors
predict_mat_remain_final <- predict_mat_remain[
  , which(predictors %in% predictors_orig), drop = FALSE
]

# TODO: can functionalise this recursively
# calculate_sigma <- function(mles, predict_mat) {
#   
# }

# sigma_t <- exp(
#   mles[1] + 
#     mles[2] * predict_mat_remain_final[, 1]
#     # mles[4] * predict_mat_remain[, 2]
#     # mles[5] * predict_mat_remain[, 3] + 
#     # mles[6] * predict_mat_remain[, 4] + 
#     # mles[7] * predict_mat_remain[, 5] 
# )

sigma_t <- exp(
  mles[1] +
    mles[2] * predict_mat_remain_final[, 1] +
    mles[3] * predict_mat_remain_final[, 2] +
    mles[4] * predict_mat_remain_final[, 3] +
    mles[5] * predict_mat_remain_final[, 4]
)

# make predictions
set.seed(123)
# TODO: Try with all data, not just remaining!
predictions <- rgpd(
  # nrow(dat_remain), shape = last(best_fit$mle), scale = sigma_t
  nrow(dat_remain), shape = last(mles), scale = sigma_t, mu = threshold #  + 0.25
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

# plot cumulatively
dat_remain_pred %>% 
  select(date, abs_diff, predictions) %>% 
  # dplyr::filter(abs_diff < 90) %>% 
  pivot_longer(abs_diff:predictions) %>% 
  group_by(name) %>% 
  mutate(value = cumsum(value)) %>% 
  ggplot(aes(x = date, y = value, colour = name)) + 
  ggtitle(paste0("Threshold = ", round(threshold, 3))) + 
  # geom_line() + 
  geom_point() + 
  theme_bw()

# Doesn't look too bad! Systematically lower prediction
# However, seems fixable! Captures temporal trend quite well!
# Seasonal temporal effect would be useful!
# Also, add uncertainty!!
  
# R^2 only maxes sense for constant variance, 
# for GPD variance varies quite a bit!!
# Use WAIC, AIC, etc, likelihood ratio test
# Q-Q plot
# Assessing predictions: Use return levels, for multivariate we need to 
# specify covariates and sample from them using MCMC (look into this!)
# But I'm only looking at one location with non-stationary extremes, so this 
# should be fine http://www.mas.ncl.ac.uk/~nlf8/teaching/mas8391/background/chapter5.pdf
# Also sample from parameter values
