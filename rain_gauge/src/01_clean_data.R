#### Load, aggregate and explore Rain Gauge, Satellite and Sonde Data ####

# TODO:
# Rain gauge:
# Load data (done)
# Identify what each column refers to!
# Take averages and sums for atmospheric conditions by day (done)
# Plot correlation between rain and other variables (done)
# Do some analysis on possible seasonal patterns in rain (still todo)
# Save (done)

# Radar:
# Take lon and lat nearest to gauge (-97.6242, 36.6656) (done)
# - On this, should I do some kringing instead?
# Take sums for rain by day (done)
# Save (done)

# Sonde: 
# Take lon and lat nearest to gauge
# Take averages and sums for atmospheric conditions by day and height
# Look at correlations with precipitation, at what height is it the most?
# Could plot atmospheric conditions vs heights, smooth across different times
# Look into what Anakin did for weather heights (done)
# Save


#### libs ####

library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
# for loading xarray based .nc data
library(ncdf4)  
library(ncdf4.helpers)
library(parallel)
# library(arrow)
library(ggcorrplot)
library(data.table)

source("src/functions.R")

#### Metadata ####

cores <- min(5, detectCores() - 1)

dir <- "rain_gauge/data/"

# already ran for 2022/2023 for gauge
gauge_files <- paste0(
  # dir, "/2022_2023/rain_gauge_daily_2022_2023.csv.gz"
  dir, "/rain_gauge_daily_full.csv.gz"
)

# gauge_files <- list.files(
#   file.path(dir, "arm_precip_data/"), pattern = ".nc", full.names = TRUE
# )

# gauge_files <- list.files(
#   file.path(dir, "2022_2023/arm_precip_data_2022_2023"), pattern = ".nc", full.names = TRUE
# )
sat_files <- list.files(
  # file.path(dir, "nexrad_reduced_sgp_2022_2023/"), 
  file.path(dir, "nexrad_reduced_sgp/"), 
  pattern = ".nc", 
  full.names = TRUE
)
# one sat file is not working!
# sat_files <- sat_files[basename(sat_files) != "st4_conus.2022070123.01h.nc"]

# sat_files_missing <- list.files(
#   file.path(dir, "nexrad_reduced_sgp_missing_2022_2023/"), 
#   pattern = ".nc", 
#   full.names = TRUE
# )

# sonde_files <- list.files(
#   file.path(dir, "arm_sonde_data/"), pattern = ".nc", full.names = TRUE
# )
sonde_files <- c(
  paste0(dir, "/arm_sonde_data/arm_sonde_binned.nc"),
  paste0(dir, "2022_2023/arm_sonde_data_2022_2023/arm_sonde_binned.nc")
)

#### Rain Gauge ####

# load rain gauge data
if (length(gauge_files) > 1) {
  gauge <- bind_rows(mclapply(gauge_files, load_gauge_netcdf, mc.cores = cores))
  
  # gauge_unlist
  
  # write to file
  # readr::write_csv(gauge, "rain_gauge/data/rain_gauge_2022_2023.csv.gz")
  # readr::write_csv(
  #   # x = data.frame(gauge), 
  #   x = data.frame(gauge, stringsAsFactors = FALSE),
  #   # file = paste0(dir, "rain_gauge_2022_2023.csv.gz")
  #   file = "rain_gauge/data/rain_gauge_2022_2023.csv.gz"
  # )
  
  # variables to take the total of
  sum_vars <- c(
    # "wxt_precip_rate_mean",
    # "wxt_cumul_precip",
    "tbrg_precip_total"
  )
  
  gauge_daily <- gauge %>% 
    mutate(date = as_date(substr(time, 0, 10))) %>% 
    group_by(date, lat, lon, alt) %>% 
    summarise(
      # take mean for certain variables
      across(all_of(mean_vars), ~ mean(.x, na.rm = TRUE)),
      across(all_of(mean_vars), ~ min(.x, na.rm = TRUE), .names = "{col}_min"),
      across(all_of(mean_vars), ~ min(.x, na.rm = TRUE), .names = "{col}_max"),
      # take sum for certain variables
      across(all_of(sum_vars), ~ sum(.x, na.rm = TRUE)), 
      .groups = "drop"
    )
 
  
  # Variables etc:
  # time, lon, lat, alt: Obvious
  # qc_x - 
  # rh - 
  # arith_mean vs vec_mean - 
  # x_std - before taking mean, can safely remove
  # wct_x - 
  # tbrg_x - 
  # logger_volt - 
  
  # understand what variables are 
  gauge %>% 
    # just take one day, to see if variables are means or cumulative
    filter(grepl("2022-07-29", as.character(time))) %>% 
    relocate(lat, lon, alt, .after = "time") %>% 
    pivot_longer(
      temp_mean:qc_logger_temp, values_to = "value", names_to = "name"
    ) %>% 
    mutate(
      name = factor(
        name, 
        levels = names(gauge[!names(gauge) %in% c("time", "alt", "lon", "lat")])
      )
    ) %>% 
    ggplot(aes(x = time, y = value)) +
    geom_line() + 
    facet_wrap(~ name, scales = "free") + 
    theme(
      axis.text.x = element_text(angle = -45, vjust = 0.5)
    )
  
  # variables to take the mean of (also take max and min)
  mean_vars <- c(
    "temp_mean",
    "rh_mean",
    "atmos_pressure",
    "wspd_arith_mean",
    "wspd_vec_mean",
    "wdir_vec_mean",
    "logger_volt", 
    "logger_temp"
  )
  
  # variables to take the total of
  sum_vars <- c(
    # "wxt_precip_rate_mean",
    # "wxt_cumul_precip",
    "tbrg_precip_total"
  )
  
  # take mean/sum of atmospheric quantities across all days, where appropriate
  # should I actually take max and min as well?
  gauge_daily <- gauge %>% 
    mutate(date = as_date(substr(time, 0, 10))) %>% 
    group_by(date, lat, lon, alt) %>% 
    summarise(
      # take mean for certain variables
      across(all_of(mean_vars), ~ mean(.x, na.rm = TRUE)),
      across(all_of(mean_vars), ~ min(.x, na.rm = TRUE), .names = "{col}_min"),
      across(all_of(mean_vars), ~ min(.x, na.rm = TRUE), .names = "{col}_max"),
      # take sum for certain variables
      across(all_of(sum_vars), ~ sum(.x, na.rm = TRUE)), 
      .groups = "drop"
    )
  
  # Now plot again
  gauge_daily %>% 
    pivot_longer(
      temp_mean:tbrg_precip_total, values_to = "value", names_to = "name"
    ) %>% 
    mutate(
      name = factor(
        name, 
        levels = names(gauge_daily)[
          !names(gauge_daily) %in% c("time", "alt", "lon", "lat")
        ]
      )
    ) %>% 
    ggplot(aes(x = date, y = value)) +
    geom_line() + 
    facet_wrap(~ name, scales = "free") + 
    theme(
      axis.text.x = element_text(angle = -45, vjust = 0.5)
    )
  
  # look at correlations between variables (particularly with rainfall)
  gauge_daily %>% 
    select(-c(lon, lat, alt)) %>% 
    mutate(date = as.numeric(date)) %>% 
    cor() %>% 
    ggcorrplot(
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
  # precipitation (and all other variables!) seems to be very weekly correlated
  # with other variables/possible predictors
  # May belie seasonal patterns here though
  rm(gauge); gc()
} else {
  gauge_daily <- readr::read_csv(gauge_files)
}

# take lon and lat vals to filter satellite data by
# gauge_lon <- unique(gauge_daily$lon)
gauge_lon <- -97.6
# gauge_lat <- unique(gauge_daily$lat)
gauge_lat <- 36.7
gauge_daily$lon <- gauge_lon
gauge_daily$lat <- gauge_lat

# save
# readr::write_csv(gauge_daily, paste0(dir, "rain_gauge_daily.csv.gz"))
readr::write_csv(gauge_daily, paste0(dir, "rain_gauge_daily_full.csv.gz"))

rm(gauge_daily); gc()


#### Satellite ####

# load satellite data
# sat_files_all <- c(sat_files, sat_files_missing)
sat_files_all <- sat_files
# split into 100, waaay too many files
splits <- rep(
  seq_len(100), 
  length.out = length(sat_files_all), 
  each = ceiling(length(sat_files_all) / 100)
)

# sat <- lapply(1:100, \(i) {
sat <- mclapply(1:100, \(i) {
  # print(paste0("Completed: ", x / 100 * 100, "%"))
  system(sprintf('echo "\n%s\n"', paste0(i, "% completed", collapse = "")))
  gc()
  # files to load
  sat_files_temp <- sat_files_all[splits == i]
  # Take point closest to in euclidean distance to rain gauge
  bind_rows(lapply(sat_files_temp, load_sat_netcdf)) %>% 
    mutate(
      dist_from_gauge = sqrt((lon - gauge_lon) ^ 2 + (lat - gauge_lat) ^ 2)
    ) %>% 
    filter(
      dist_from_gauge == min(dist_from_gauge)
    ) %>% 
    select(-dist_from_gauge) %>% 
    # also take cumulative precip across all days
    mutate(date = as_date(substr(time, 0, 10))) %>% 
    group_by(lon, lat, date) %>% 
    summarise(precip = sum(precip, na.rm = TRUE), .groups = "drop")
},  mc.cores = cores)

# aggregate again, for dates across multiple files
sat <- bind_rows(sat) %>% 
  group_by(lon, lat, date) %>% 
  summarise(precip = sum(precip, na.rm = TRUE), .groups = "drop")

# plot
sat %>% 
  ggplot(aes(x = date, y = precip)) + 
  geom_line()

# sat_22_23 <- readr::read_csv(
#   "rain_gauge/data/2022_2023/sat_dat_daily_2022_2023.csv.gz"
# )
# sat <- sat %>% 
#   bind_rows(sat_22_23) %>% 
#   arrange(date)

# save and remove
readr::write_csv(sat, paste0(dir, "sat_dat_daily.csv.gz"))

rm(sat); gc()

#### Sonde ####

# TODO: Make this like sat data, much quicker!
if (length(sonde_files) > 2) {
  sonde <- bind_rows(mclapply(seq_along(sonde_files), \(i) {
  # sonde <- bind_rows(mclapply(seq_along(sonde_files[1:100]), \(i) {
    gc()
    # give progress update
    system(sprintf(
      'echo "\n%s\n"', 
      paste0(
        round(100 * (i / length(sonde_files)), 2), "% completed", collapse = ""
      )
    ))
    load_sonde_netcdf(sonde_files[i]) %>% 
      # average over days and height
      mutate(date = as_date(substr(time, 0, 10))) %>% 
      group_by(height, date) %>% 
      summarise(across(everything(), ~ mean(.x)), .groups = "drop")
  }, mc.cores = cores))
} else {
  # sonde <- load_sonde_netcdf(sonde_files, height_binned = TRUE)
  sonde <- bind_rows(lapply(
    sonde_files, load_sonde_netcdf, height_binned = TRUE
  ))
  
  # temp: do 2022/2023 sonde data
  sonde <- sonde[[2]]
}
gc()


# average again over days and height, in case different files have same times
# sonde_daily <- sonde %>% 
#   mutate(date = as_date(substr(time, 0, 10))) %>% 
#   select(-time) %>% 
#   group_by(height, date) %>% 
#   summarise(across(everything(), ~ mean(.x)), .groups = "drop") %>% 
#   mutate(height = as.numeric(height))

# too slow! Use data.table
# sonde_daily <- setDT(sonde)
# sonde_daily <- sonde_daily[, date := as.Date(substr(time, 1, 10))][, 
#   time := NULL
# ][, 
#   lapply(.SD, mean), by = .(height, date)
# ][, 
#   height := as.numeric(height)
# ]

sonde_dt <- setDT(sonde)
# split in three, too computationally expensive otherwise
# TODO: Split by year instead!?
n_splits <- 10
sonde_list <- split(
  sonde_dt, 
  rep(
    seq_len(n_splits), 
    length.out = nrow(sonde), 
    each = ceiling(nrow(sonde) / n_splits)
  )
)

rm(sonde, sonde_dt); gc()

i <- 0
sonde_daily <- lapply(sonde_list, \(x) {
  i <<- i + 1
  print(paste0(i / n_splits * 100, "% finished"))
  gc()
  x[, date := as.Date(substr(time, 1, 10))][, 
    time := NULL
  ][, 
    lapply(.SD, mean), by = .(height, date)
  ][, 
    height := as.numeric(height)
  ]
})

# Aggregate once more
sonde_daily <- rbindlist(sonde_daily)
sonde_daily <- sonde_daily[,
  lapply(.SD, mean), by = .(height, date)
]

rm(sonde_list); gc()

# plot against height and against time
sonde_daily_long <- sonde_daily %>% 
  # pivot_longer(bar_pres:wspd)
  pivot_longer(rh:wspd)

# Seems to be pretty constant over time, but w direction increased in 2022/23
p_time <- sonde_daily_long %>% 
  ggplot(aes(x = date, y = value)) + 
  geom_point() + 
  geom_smooth(method = "gam") + 
  facet_wrap(~ name)

# Seems to pretty constant as well, doesn't really justify splitting by height
p_height <- sonde_daily_long %>% 
  ggplot(aes(x = height, y = value)) + 
  geom_point() + 
  geom_smooth(method = "gam") + 
  facet_wrap(~ name)

# sonde_22_23 <- readr::read_csv(
#   "rain_gauge/data/2022_2023/sonde_daily_2022_2023.csv.gz"
# )
# sonde_daily <- sonde_daily %>% 
#   bind_rows(sonde_22_23) %>% 
#   arrange(date)

sonde_daily_no_height <- sonde_daily[,
  lapply(.SD, mean), by = .(date)
]

# save
# readr::write_csv(sonde_daily, paste0(dir, "sonde_daily.csv.gz"))
readr::write_csv(sonde_daily, paste0(dir, "sonde_daily_22_23.csv.gz"))
readr::write_csv(
#   sonde_daily_no_height, paste0(dir, "sonde_daily_no_height.csv.gz")
  sonde_daily_no_height, paste0(dir, "sonde_daily_no_height_22_23.csv.gz")
)

rm(sonde, sonde_daily, sonde_daily_long); gc()
