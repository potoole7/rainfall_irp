#### Initial Exploration of Precipitation Data ####

# TODO: Look at relative difference as well!! Just exploratory

# - Data only goes up to April 2023, is that correct?
# - Break in data in Nov 2022 and 2023-04-01, where there is no sat precip
# - Two of those November dates coincide with quite heavy gauge rain, can't 
#   afford to be missing data!

#### libs ####

library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggcorrplot)

library(ecd)
library(ismev)

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



#### Merge datasets ####

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
  )

#### Explore Rainfall patterns ####

# missing data in gauge & satellite data
all_dates <- seq.Date(min(dat$date), max(dat$date), by = 1)
# data complete from April 2022 to April 2023, need 2017-2021 from Marion!
missing_gauge_dates <- all_dates[!all_dates %in% gauge$date]
# Missing 2019-04-04, November 2022, 2023-04-01
missing_sat_dates <- all_dates[!all_dates %in% sat$date]

print(paste0("Days of data: ", nrow(dat))) # 2191/366
print(paste0(
  "Days with no rain for either dataset: ", 
  dat %>% 
    filter(tbrg_precip_total == 0, sat_precip == 0) %>% 
    nrow()
)) # 240 (not that bad!) 

print(paste0(
  "Days with more rain for satellite: ", 
  dat %>% 
    filter(more_rain == "satellite") %>% 
    nrow()
)) # 89
print(paste0(
  "Days with rain for satellite, but not for gauge: ",
  dat %>% 
    filter(tbrg_precip_total == 0, sat_precip > 0) %>% 
    nrow()
)) # 33 => 56 days where both have rain but satellite has more, lots of info!

print(paste0(
  "Days with more rain for rain gauge: ", 
  dat %>% 
    filter(more_rain == "gauge") %>% 
    nrow()
)) # 6 (much rarer event, may mean it's easier to just take abs diff?)
print(paste0(
  "Days with rain for gauge, but not for satellite: ",
  dat %>% 
    filter(tbrg_precip_total > 0, sat_precip == 0) %>% 
    nrow()
)) # 3

# dates with NAs for diff
# date_seq <- seq(min(dat$date), max(dat$date), by = 1)
# (missing_dates <- date_seq[!date_seq %in% dat[!is.na(dat$diff), ]$date])


#### Plots ####

# plot time series
dat %>% 
  select(
    date, 
    `Gauge Precipitation` = tbrg_precip_total, 
    `Satellite Precipitation` = sat_precip
  ) %>% 
  pivot_longer(2:3) %>% 
  ggplot(aes(x = date, y = value)) + 
  geom_point() + 
  facet_wrap(~ name)
# does seem to be one instance with much more gauge rain, much larger 
# discrepancy for radar though!

dat %>% 
  mutate(
    year = substr(date, 0, 4), 
    # ignore warning, due to leap year date
    date = as_date(paste0("2022-", substr(date, 6, 10)))
  ) %>% 
  select(
    date,
    `Gauge Precipitation` = tbrg_precip_total, 
    `Satellite Precipitation` = sat_precip,
    year
  ) %>% 
  pivot_longer(2:3) %>%
  ggplot(aes(x = date, y = value, colour = name)) + 
  geom_line(alpha = 0.5, size = 1) + 
  facet_wrap(~ year, ncol = 1) +
  scale_x_date(labels = scales::date_format("%b"), date_breaks = "month") + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))
# does seem to be one instance with much more gauge rain, much larger 
# discrepancy for radar though!


dat %>% 
  filter(!is.na(diff)) %>% 
  select(
    date, `Precipitation Difference` = diff, `Absolute Difference` = abs_diff
  ) %>% 
  pivot_longer(2:3) %>% 
  ggplot(aes(x = date, y = value)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~ name)
# does seem to be one instance with much more gauge rain, much larger 
# discrepancy for radar though!

# Also plot correlation between diff and atmospheric variables
# rearrange in alphabetical order
dat[, order(colnames(dat))] %>% 
  select(-c(lon, lat, alt, more_rain)) %>% 
  relocate(date, diff, abs_diff) %>% 
  mutate(date = as.numeric(date)) %>% 
  filter(!is.na(diff)) %>% 
  cor() %>% 
  ggcorrplot(
    outline.color = "white",
    ggtheme = theme_bw(),
    tl.cex = 8, 
    legend.title = "Correlation"
  ) + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 15, hjust = 0.5)
  )
# Again, doesn't seem to be particularly correlated with much
# Only things it appears to be correlated to are things it may be correlated to:
# date, rh_mean (min, max), atmos_pressure (min, max)
# These variables are themselves correlated to others, so may be interactions
#' *Also, need to get weather data from sonde, more accurate!*

#### Test for correlations ####

# All data
dat[complete.cases(dat), ] %>% 
  select(-c(lon, lat, alt, more_rain)) %>% 
  mutate(date = as.numeric(date)) %>% 
  cor() %>% 
  ggcorrplot::ggcorrplot(
    outline.color = "white",
    # show.legend = has_legend[i],
    ggtheme = theme_bw(),
    tl.cex = 8, 
    # title = site_names[i],
    legend.title = "Correlation" # ,
    # method = "circle"
  ) + 
  # scale_x_discrete(labels = labels) + 
  # scale_y_discrete(labels = labels) + 
  # theme_bw(base_size = 12) + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 15, hjust = 0.5)
  )

# differences over some amount
dat[complete.cases(dat), ] %>% 
  filter(abs_diff > 20) %>% 
  select(-c(lon, lat, alt, more_rain)) %>% 
  mutate(date = as.numeric(date)) %>% 
  cor() %>% 
  ggcorrplot::ggcorrplot(
    outline.color = "white",
    # show.legend = has_legend[i],
    ggtheme = theme_bw(),
    tl.cex = 8, 
    # title = site_names[i],
    legend.title = "Correlation" # ,
    # method = "circle"
  ) + 
  # scale_x_discrete(labels = labels) + 
  # scale_y_discrete(labels = labels) + 
  # theme_bw(base_size = 12) + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 15, hjust = 0.5)
  )

#### Extremes Investigation ####

# Mean residual life plot
par(mfrow = c(1, 1))
evd::mrlplot(dat$abs_diff) # any threshold between > 0 and ~12 seems reasonable

# Threshold Choice Plot
par(mfrow = c(1, 2))
evd::tcplot(dat$abs_diff, tlim = c(2, 5))
evd::tcplot(dat$abs_diff, tlim = c(12, 15))

# fit simple GPD model
fit <- gpd.fit(xdat = dat$abs_diff[!is.na(dat$abs_diff)], threshold = 10)
# 17 excesses for threshold of 10, 12 for 12, hopefully more data will help!

# return level plot gives 40 exceedances in 1 year! Only observe 17 vals above 10
# in one year of data
gpd.diag(fit)





