# install.packages("terra", lib = "/projappl/project_2003061/Rpackages/")
# Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
# remotes::install_github("rspatial/terra", lib = "/projappl/project_2003061/Rpackages/")
# install.packages("rgee", lib = "/projappl/project_2003061/Rpackages/")
# remotes::install_github("r-spatial/rgee", lib = "/projappl/project_2003061/Rpackages/")
# remotes::install_github("rstudio/reticulate", lib = "/projappl/project_2003061/Rpackages/")
library(reticulate)
library(rgee, lib.loc = "/projappl/project_2003061/Rpackages/")
library(sf)
library(tidyverse)
library(lubridate)
library(foreach)
library(doParallel)
library(R.utils)
library(terra)

# rgee::ee_clean_pyenv()
# rgee::ee_install_set_pyenv(py_path = "/opt/rh/rh-python38/root/usr/bin/python")

ee_Initialize("pekkaniittynen", drive = T)

# If authentication failes, try in Puhti computing node
# module load r-env-singularity
# echo "PATH=$PATH:/appl/opt/google-cloud-sdk/bin/" >> ~/.Renviron
# start-r
# library(rgee, lib.loc = "/projappl/project_2003061/Rpackages/")
# ee_Initialize()

# Check if everything looks fine with rgee
reticulate::py_available()
ee_check()
ee_user_info()
ee_users()

# Working directory has to be specified if the script is run in a non-interactive session
# setwd("/projappl/project_2003061/repos/Sentinel_Rastigaisa")
source("scripts/FUNCTION_extract_sentinel2_NDVI_gee.R")
source("scripts/FUNCTION_check_raster.R")

# Detect how many cores are available
workers <- min(c(length(future::availableWorkers()),
                 future::availableCores()))
print(workers)

# Set a "project name" e.g., the name of the study area. Will be used in naming the files
iAREA <- "aland"

my_timeout <- 2000 # in seconds, maximum time to wait for downloading the image file from GEE
info_timeout <- 20 # in seconds, maximum time to wait for image info from GEE

mindate <- "2017-05-01"
maxdate <- Sys.Date()

# Specify the EPSG code for wanted projection
epsg <- 32634

# Sequence of dates to query the imagery
dates <- as.character(seq.Date(as.Date(mindate), as.Date(maxdate), by = "days"))
dates <- dates[month(dates) %in% 5:8] # Limit which months will be included

# Direction where imagery will be stores (in area-specific sub-directories)
base_dl_dir <- "/scratch/project_2003061/microclim_RS/sentinel2"

# Read in the study area polygon

aoi <- st_read("gis/aland_aoi_polygon.gpkg") %>%
  st_transform(crs = epsg)

# AOI polygon to GEE format
aoi_ee <- aoi %>% st_transform(crs = 4326) %>% 
  st_geometry() %>% 
  sf_as_ee()

# Create sub-directories in Google Drive and local disc where imagery will be saved if does not exists
drive_folder <- paste0("TEMP_S2_",iAREA)
area_dl_dir <- paste0(base_dl_dir,"/",iAREA)
if(!dir.exists(area_dl_dir)){
  dir.create(area_dl_dir)
}

# First round of downloads
# List already downloaded images if any
tifs <- list.files(area_dl_dir, pattern = "_S2_.*.tif$")

image_df <- tibble(file = tifs)
image_df$date <- ymd(substr(unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][5])), 1, 8))

# Remove dates with already downloaded imagery
dates2 <- dates[!dates %in% as.character(image_df$date)]

lss <- mclapply(dates2, extract_sentinel2_NDVI_gee, mc.cores = workers,
                aoi_ee = aoi_ee, epsg = epsg, drive_folder = drive_folder, 
                area_dl_dir = area_dl_dir, name = iAREA)

ee_clean_container(name = drive_folder, type = "drive", quiet = FALSE)

# Second round of downloads
# List downloaded images

tifs <- list.files(area_dl_dir, pattern = "_S2.*.tif$")

image_df <- tibble(file = tifs)
image_df$date <- ymd(substr(unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][5])), 1, 8))

dates2 <- dates[!dates %in% as.character(image_df$date)]
dates2 <- sample(dates2, length(dates2))

lss <- mclapply(dates2, extract_sentinel2_NDVI_gee, mc.cores = workers,
                aoi_ee = aoi_ee, epsg = epsg, drive_folder = drive_folder, 
                area_dl_dir = area_dl_dir, name = aoi$name)

ee_clean_container(name = drive_folder, type = "drive", quiet = FALSE)

# List downloaded images
tifs <- list.files(area_dl_dir, pattern = "_S2.*.tif$")

image_df <- tibble(file = tifs)
image_df$area <- unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][1]))
image_df$satid <- unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][2]))
image_df$proclevel <- unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][3]))
image_df$datatype <- unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][4]))
image_df$tile <- gsub(".tif","",unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][7])))
image_df$date <- ymd(substr(unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][5])), 1, 8))
image_df$time <- substr(unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][5])), 10, 15)

image_df <- image_df %>% arrange(date)

# Check if all rasters are working. Remove corrupted files.
img_remove <- unlist(mclapply(image_df$file, check_raster, image_dir = area_dl_dir, mc.cores = workers))

if(length(img_remove) > 0){
  print(paste0(length(img_remove), " raster(s) not functional. Try download again!!"))
  dates3 <- as.character(ymd(substr(unlist(lapply(image_df$file, function(x) str_split(x, "_")[[1]][5])), 1, 8)))
  lss <- mclapply(dates3, extract_sentinel2_NDVI_gee, mc.cores = workers,
                  aoi_ee = aoi_ee, epsg = epsg, drive_folder = drive_folder, 
                  area_dl_dir = area_dl_dir, name = aoi$name)
}

img_remove <- unlist(mclapply(image_df$file, check_raster, image_dir = area_dl_dir, mc.cores = workers))

if(length(img_remove) > 0){
  print(paste0(length(img_remove), " raster(s) not functional. REMOVED!!"))
  unlink(paste0(area_dl_dir,"/",img_remove))
}

image_df <- image_df %>% 
  filter(!file %in% img_remove)

image_df %>% write_csv(paste0(area_dl_dir,"/image_info.csv"))

unlink(list.files(tempdir(), full.names = T))
ee_clean_container(name = drive_folder, type = "drive", quiet = FALSE)
