#*******************************************************************************
# Project: roboBayes pipeline
# Script Purpose: plotting harmonic models from roboBayes results
# Date: 2022-12-22
# Author: Jenna Abrahamson
#*******************************************************************************
library(terra)
library(fields)
library(inlmisc)
library(optparse)
library(RJSONIO)
library(data.table)
source("../functions/utils.R")
source("../functions/post_extract_quantities.R")
source("../functions/post_filters.R")


# #############################################################################
# Load all data
# #############################################################################

# read in og data
og_data <- "./ACC.Rdata"
load(og_data)

# read in results
results_path <- "./US_R004_11SKD57_ndvi_swir2_high_low_soil_qa.Rdata"
load(results_path)

tmp_s_path <- "./tmp.tif"
tmp_s <- rast(tmp_s_path)$tmp_1
values(tmp_s) <- NA
values(tmp_s) <- c(t(matrix((1:length(results)),nrow(tmp_s),ncol(tmp_s))))
writeRaster(tmp_s, "./US_R004_57_Index.tif")


nsample <- 10000
season <- 1
bands <- c("ndvi", "swir2", "high", "low", "vege", "soil", "qa", "nir", "green")
model_type <- "same"
vscale <- 1
possible_signals <- c("blue", "green", "red", "nir", "swir1", "swir2",
    "ndvi", "high", "low", "vege", "soil", "qa", "brighness")



##### We can probably combine these steps
# necessary to format data correctly
data_cube <- input_data_check(data_cube, dates, index)
# scale data to be between 0 and 1
data_cube <- signal_scaling(data_cube, possible_signals)
# Elimanate data that does not obey the data bounds
data_cube <- domain_filter(data_cube, possible_signals)
# Eliminate data contaminated by clouds, aerosols, etc.
data_cube <- quality_filter(data_cube)
# Eliminate data with water (optional)
data_cube <- water_filter(data_cube)
# Retain only the data associated with the satellite source of interest
# data_cube <- source_filter(data_cube, sensors = sensors, sat = sat)

data_cube <- data_cube[,,unlist(bands)]
# #############################################################################
# Post-Processing and Plotting
# #############################################################################

#### call whatever functions we want to get summarized data
all_cps <- get_all_cps(results, dates)
cps_first <- get_first_cps(results, dates)
cps_latency <- get_cps_min_latency(results, dates)

# Make changepoint raster
cps_first_rast <- rast(tmp_s)
values(cps_first_rast) <- c(t(matrix(cps_first,nrow(tmp_s),ncol(tmp_s))))


#### Get initial polygons
polys <- to_polygons(cps_first_rast)

#### Filters

# small site filter - there should maybe be a delay in when the small site filter
# is applied so as not to lose beginnings of sites
polys <- (small_site_filter(polys, size_thresh = 9000))

# make end date the maximum
polys_filtered$end_date <- max(dates)
polys_filtered <- project(polys_filtered, "epsg:4326")


######################################
# Plotting
######################################

plot_all <- function(pix, data_cube, results, dates){
    pix <- pix
    ndvi_pix <- as.data.frame(data_cube[pix,,"ndvi"])
    swir_pix <- as.data.frame(data_cube[pix,,"swir2"])
    soil_pix <- as.data.frame(data_cube[pix,,"soil"])
    high_pix <- as.data.frame(data_cube[pix,,"high"])
    low_pix <- as.data.frame(data_cube[pix,,"low"])
    vege_pix <- as.data.frame(data_cube[pix,,"vege"])

    par(mfrow=c(6,1))
    plot(dates, ndvi_pix[,1], ylab = "NDVI", xlab="Dates", main="NDVI", pch=19)
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2.5)
        lines(c(1:52), NDVI_60512[,2], lwd=2, add=T)
    }
    title(paste0("Pixel: ", pix), cex=2, line=3)

    plot(dates, swir_pix[,1], ylab = "SWIR", xlab="Dates", main="SWIR2", pch=19)
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2.5)
        
    }
    plot(dates, soil_pix[,1], ylab = "Soil", xlab="Dates", main="SOIL", pch=19)
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2.5)
    }

    plot(dates, high_pix[,1], ylab = "High Albedo", xlab="Dates", main="HIGH ALBEDO", pch=19)
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2.5)
    }

    plot(dates, low_pix[,1], ylab = "Low Albedo", xlab="Dates", main="LOW ALBEDO", pch=19)
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2.5)
    }

    plot(dates, vege_pix[,1], ylab = "Vegetation", xlab="Dates", main="VEGETATION", pch=19)
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2.5)
    }
}

plot_all(20381, data_cube, results, dates)

# Arguments for plotting harmonic model
before_sin_val <- results[[60512]]$mods[[1]]$B[2,1]
before_cos_val <- results[[60512]]$mods[[1]]$B[3,1]
before_mean_val <- results[[60512]]$mods[[1]]$B[1,1]

# Function for harmonic models
plot_har <- function(sin_val, cos_val, mean_val, dates) {
  Amp <- CalcAmp(sin_val, cos_val)
  x <- seq(52, (length(dates)-1), by = 1)
  sin_curve <- Amp * sin(x)
  cos_curve <- Amp * cos(x)

  har_curve <- mean_val + sin_val*sin(2*pi*x/365) + cos_val*cos(2*pi*x/365) 
  plotting_vals <- cbind(x, har_curve)
  return(plotting_vals)
}

dates_seq <- seq(1, 52, by=1)
plot(dates_seq, NDVI_60512[,2], col = "#03045E", type="l", lwd=2)

NDVI_60512 <- plot_har(before_sin_val, before_cos_val, before_mean_val)
