#### Download and preprocess Rain Gauge NetCDFs ####

#### Libs ####

library(ncdf4)
library(raster) 
library(ggplot2)
library(parallel)
library(lubridate)

# TODO: Source other files

#### Download Files ####

# load sat data
sat <- readr::read_csv("rain_gauge/data/sat_dat_daily.csv.gz")

# Specify the URL of the file you want to download
# url <- "https://archive.arm.gov/orders/fileServer/orders/otoolep1/246248/sgpwbpluvio2C1.a1.20240312.000000.nc"
url <- "https://archive.arm.gov/orders/fileServer/orders/otoolep1/246248/sgpwbpluvio2C1.a1."

# dates with files for radar data
dates <- seq.Date(min(sat$date), max(sat$date), by = "day")
dates <- stringr::str_remove_all(dates, "-")

urls <- paste0(url, dates, ".000000.nc")

# temp
# urls <- urls[length(urls)]

# Specify the file name and location where you want to save the file on your computer
# file_name <- "my_data.nc"
file_names <- basename(urls)
# file_path <- "/path/to/save/folder/"
file_path <- "~/temp/" # TODO: Change this
files <- paste(file_path, file_names, sep = "")

# Call the download.file() function, passing in the URL and file name/location as arguments
# lapply(seq_along(urls), \(i) { 
#   download.file(urls[i], files[i], mode = "wb")
# })
mc.cores <- detectCores() - 1
mclapply(seq_along(urls), \(i) {
  download.file(urls[i], files[i], mode = "wb")
})

# check that all files have been downloaded
files_success <- paste0(file_path, list.files(file_path))
files_remaining <- files[!files %in% files_success]
urls_remaining <- urls[!files %in% files_success]

# try again for remaining files
download.file(urls_remaining[1], files[1], mode = "wb")


#### Preprocess ####

nc_data <- nc_open(files[1])
