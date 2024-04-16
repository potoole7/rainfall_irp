#### Preprocess Rain Gauge NetCDFs ####

# Have five years of gauge data, want to compare to original 1 year of data 
# sent over by Marion

#### Libs ####

library(ncdf4)
library(raster) 
library(dplyr)
library(ggplot2)
library(parallel)
library(lubridate)

# TODO: Source other files
source("src/functions.R")

#### Metadata ####

dir <- "rain_gauge/data/rain_gauge_full/"
cores <- detectCores() - 1

#### Load Data ####

files <- list.files(dir, full.names = TRUE)
gauge <- bind_rows(mclapply(files, load_gauge_netcdf, mc.cores = cores))

# also load one netcdf, to understand what columns represent
nc_gauge <- nc_open(files[1])

var_descriptions <- lapply(nc_gauge$var, \(x) {
  x$longname
}) %>% 
  tibble::enframe() %>% 
  mutate(value = unlist(value))

# plot 
gauge %>% 
  select(time, accum_total_nrt) %>% 
  ggplot(aes(x = time, y = accum_total_nrt)) +
  geom_line()
  
#### Take Daily Averages ####

mean_vars <- c(
  "accum_rtnrt", "accum_nrt", "accum_total_nrt", "bucket_rt", "bucket_nrt"
)

gauge_daily <- gauge %>% 
  mutate(date = as_date(substr(time, 0, 10))) %>% 
  group_by(date, lat, lon, alt) %>% 
  summarise(
    # take mean for certain variables
    across(all_of(mean_vars), ~ mean(.x, na.rm = TRUE)),
    # across(all_of(mean_vars), ~ min(.x, na.rm = TRUE), .names = "{col}_min"),
    # across(all_of(mean_vars), ~ min(.x, na.rm = TRUE), .names = "{col}_max"),
    # take sum for certain variables
    # across(all_of(sum_vars), ~ sum(.x, na.rm = TRUE)), 
    .groups = "drop"
  )

#### Compute daily accum_total_nrt ####

# TODO: Look into doing so 


#### Compare with previous data ####

# load original data
gauge_orig_daily <- readr::read_csv("rain_gauge/data/rain_gauge_daily.csv.gz")

# join and plot
gauge_daily %>% 
  left_join(select(gauge_orig_daily, date, tbrg_precip_total), by = "date") %>% 
  # only keep dates for which we have tbrg_precip_total
  dplyr::filter(date %in% gauge_orig_daily$date) %>% 
  tidyr::pivot_longer(accum_rtnrt:tbrg_precip_total) %>% 
  dplyr::filter(
    # name %in% c("accum_nrt", "bucket_rt", "bucket_nrt", "tbrg_precip_total")
    # name %in% c("accum_nrt", "bucket_nrt", "tbrg_precip_total")
    # name %in% c("accum_nrt", "bucket_nrt", "tbrg_precip_total")
    name %in% c("tbrg_precip_total", "accum_total_nrt")
  ) %>% 
  ggplot(aes(x = date, y = value, colour = name)) + 
  geom_line(alpha = 0.5, size = 2) + 
  # scale_y_continuous(limits = c(0, 50)) + 
  # scale_x_date(date)
  NULL

#### Compare daily observations ####

gauge_files <- list.files(
  "rain_gauge/data/2022_2023/arm_precip_data_2022_2023", pattern = ".nc", full.names = TRUE
)

gauge_orig <- bind_rows(mclapply(gauge_files, load_gauge_netcdf, mc.cores = cores))

# plot daily original readings vs new readings
gauge_orig <- gauge_orig %>% 
  select(time, tbrg_precip_total)

gauge_compare <- gauge %>% 
  dplyr::filter(time %in% gauge_orig$time) %>% 
  select(time, accum_total_nrt) %>% 
  left_join(gauge_orig, by = "time") %>% 
  tidyr::pivot_longer(accum_total_nrt:tbrg_precip_total)

gauge_compare %>% 
  ggplot(aes(x = time, y = value, colour = name)) +
  geom_line()

#### Parsing New Rain Data ####

# Idea: rain data accumulates over time, then is reset when gauge is emptied
# Want to avoid these resets
# Take differences over different time points, ignore negative differences?

# Problem: Dont' have full times for each date? Explore!
# Just ignore!!
# times <- seq(from = min(gauge$time), to = max(gauge$time), by = 60)
# missing_times <- times[!times %in% gauge$time]
# # any missing days? Seems to be missing times for many dates!
# unique(as_date(missing_times))
# 
# test <- gauge %>% 
#   mutate(date = as_date(time)) %>% 
#   dplyr::filter(date == as_date("2022-02-04"))
# test_times <- seq(from = min(test$time), to = max(test$time), by = 60)
# sort(test_times[!test_times %in% test$time])


# parse cumulative rain from periodically "emptied" rain gauge
cumul_gauge <- gauge %>% 
  select(time, accum_total_nrt) %>% 
  # fill NAs with previous values
  tidyr::fill(accum_total_nrt, .direction = "down") %>% 
  mutate(
    # remove NAs (DON'T DO THIS, MAKES DIFFERENCES MUCH LARGER!)
    # accum_total_nrt = ifelse(is.na(accum_total_nrt), 0, accum_total_nrt),
    # take difference between rain at time points
    lag_accum_total_nrt = lag(accum_total_nrt, default = 0),
    # lag_accum_total_nrt = ifelse(
    #   is.na(lag_accum_total_nrt), 0, lag_accum_total_nrt
    # ),
    diff_accum = accum_total_nrt - lag_accum_total_nrt, 
    # ignore negative differences, indicating the gauge is emptied
    diff_accum = ifelse(diff_accum < 0, 0, diff_accum),
    # take cumulative differences, giving cumulative rain 
    cum_accum_total_nrt = cumsum(diff_accum),
  ) %>% 
  select(-c(lag_accum_total_nrt, diff_accum))

# plot to check things have been done correctly
# Bug: Seems to work fine for before 2020, but breaks there, investigate!
cumul_gauge %>%  
  slice(1:3e5) %>% 
  # slice(2e5:4e5) %>% 
  tidyr::pivot_longer(accum_total_nrt:cum_accum_total_nrt) %>% 
  ggplot(aes(x = time, y = value, colour = name)) +
  geom_line()

# convert from cumulative rain to daily
gauge_daily <- cumul_gauge %>% 
  mutate(
    # convert to minute
    lag_rain = lag(cum_accum_total_nrt, default = 0), 
    rain = cum_accum_total_nrt - lag_rain, 
    date = as_date(time)
  ) %>% 
  select(-lag_rain) %>% 
  # convert to daily
  group_by(date) %>% 
  summarise(rain = sum(rain), .groups = "drop")
  
# plot
gauge_daily %>% 
  ggplot(aes(x = date, y = rain)) +
  geom_point()

#### Compare to old gauge data and satellite ####

# original gauge data
gauge_orig_daily <- readr::read_csv(
  "rain_gauge/data/2022_2023/rain_gauge_daily_2022_2023.csv.gz"
)

# RADAR data
satellite_dat <- readr::read_csv("rain_gauge/data/sat_dat_daily.csv.gz")

# compare again to other gauge readings (might be easier to look cumulatively?)
gauge_compare <- gauge_daily %>% 
  # dplyr::filter(date %in% gauge_orig_daily$date) %>% 
  select(date, rain) %>% 
  left_join(
    select(gauge_orig_daily, date, old_gauge = tbrg_precip_total), 
    by = "date"
  ) %>% 
  left_join(select(satellite_dat, date, radar_precip = precip)) %>% 
  rename(new_gauge = rain) %>% 
  # don't use extra dates in new gauge data!
  dplyr::filter(date %in% satellite_dat$date)

dot_plot <- function(dat, filter_zeros = TRUE) {
  x <- dat %>% 
    # tidyr::pivot_longer(rain:tbrg_precip_total) %>% 
    # filter(date %in% gauge_orig_daily$date) %>% 
    tidyr::pivot_longer(new_gauge:radar_precip)
  
  if (filter_zeros) x <- dplyr::filter(x, !value == 0)
  x %>% 
    ggplot(aes(x = date, y = value, colour = name)) + 
    geom_point()
}

dot_plot(gauge_compare)

dot_plot(dplyr::filter(gauge_compare, date %in% gauge_orig_daily$date))

# Interestingly identifies outliers for new gauge data; look further
gauge_compare %>% 
  filter(new_gauge > 90) %>% 
  mutate(date = 1:n()) %>% 
  dot_plot(filter_zeros = FALSE)
# dot_plot(dplyr::filter(gauge_compare, new_gauge > 90), filter_zeros = FALSE)
# 2018-05-01:2018-05-02: Gauge records very high rainfall not in radar data!

# cumulative
cum_line_plot <- function(dat) {
  dat %>% 
    tidyr::pivot_longer(new_gauge:radar_precip) %>% 
    # filter out possible outlier?
    # filter(value < 60) %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    group_by(name) %>% 
    mutate(value = cumsum(value)) %>% 
    ungroup() %>% 
    mutate(value = ifelse(value == 0, NA, value)) %>% 
    ggplot(aes(x = date, y = value, colour = name)) + 
    geom_line(size = 2, alpha = 0.7)
}

cum_line_plot(gauge_compare)

cum_line_plot(dplyr::filter(gauge_compare, date %in% gauge_orig_daily$date))

#### Temp: Use 22-23 weather covariates from gauge for all years ####

# take rain variables for 22-23
rain_covars <- gauge_orig_daily %>% 
  select(date, temp_mean, rh_mean, atmos_pressure, wspd_vec_mean) %>% 
  # take month and day from date
  mutate(month_day = substr(date, 6, 10))

# Use weather data for this year for all years 
years <- unique(substr(gauge_daily$date, 0, 4))
rain_covars_all <- lapply(years, function(x) {
  # rain_covars$year = x
  # print(x)
  rain_covars %>% 
    mutate(year = x)
}) %>% 
  bind_rows() %>% 
  mutate(date = as_date(paste(year, month_day, sep = "-"))) %>% 
  select(-c(year, month_day))

# expand.grid(rain_covars, years)

#### Save New Data ####

# Weird error with list date, fix with below
gauge_daily %>% 
  mutate(date = as_date(as.numeric(date))) %>% 
  readr::write_csv("rain_gauge/data/rain_gauge_daily_full.csv.gz")

readr::write_csv(rain_covars_all, "rain_gauge/data/rain_covars.csv.gz")
