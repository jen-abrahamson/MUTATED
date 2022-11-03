#*******************************************************************************
# Project: MUTATED - roboBayes
# Script Purpose: Investigate KR_R002 False Positives/Missed Change Dates
# Date: 2022-07-14
# Author: Jenna Abrahamson
#*******************************************************************************
# Load necessary functions/libraries, set path to functions scripts 
source("./robobayes/functions/data_prep.R")
source("./functions/priors_mv_reg.R")
source("./functions/scaling_mv_reg.R")
source("./functions/utils.R")
source("./functions/spike_filter.R")
source("./functions/yaml_utils.R")
library(foreach)
library(future.apply)
library(optparse)
library(nlme)
library(parallel)
library("RJSONIO")
plan(multisession)


# #############################################################################
option_list <- list(
  make_option(c("-p", "--params"),
    type = "character",
    default = "../../../cirrus/params.json",
    help = "Path to param json.",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# set path to parameter file on your machine
params <- fromJSON("./cirrus/KR_R002_T7.json")$ROBOBAYES
data_dir <- params$DATA_DIR
model_params <- params$MODEL
spike_type <- params$SPIKE
spike_params <- params$SPIKEPARAMS
sensors <- params$SENSORS
nsample <- params$NSAMPLE
change_date <- params$CHANGE_DATE
season <- params$SEASON
bands <- params$BANDS
model_type <- params$MODEL_TYPE
npix <- params$NPIX
vscale <- params$VSCALE
possible_signals <- params$POSSIBLE_SIGNALS
out_tag <- params$OUT_TAG
region <- params$REGION


# Load data
load("./ACC.Rdata")
# Load results
load("./ACC.Rdata")

plots_save <- "./Eval3Plots"
# Find number of pixels
if (is.null(npix)) {
  npix <- nrow(data_cube[, , "qa"])
}


# #############################################################################

if (is.null(change_date)) {
  change_date <- max(dates)
}

##### Do preprocessing of raw data
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
# Plotting functions
# #############################################################################

# Takes in a pixel number, data_cube, results, and dates
# Plots all time series with change points and
# x's over removed spike spoints

plot_all <- function(pix, data_cube, results, dates, plot_directory){
    pix <- pix
    ndvi_pix <- as.data.frame(data_cube[pix,,"ndvi"])
    swir_pix <- as.data.frame(data_cube[pix,,"swir2"])
    high_pix <- as.data.frame(data_cube[pix,,"high"])
    low_pix <- as.data.frame(data_cube[pix,,"low"])
    soil_pix <- as.data.frame(data_cube[pix,,"soil"])
    vege_pix <- as.data.frame(data_cube[pix,,"vege"])

    pdf(file = paste0(plot_directory, paste("/KR_R002-", pix, ".pdf", sep="")), width=9, height=16)
    par(mfrow=c(6,1))
    plot(dates, ndvi_pix[,1], ylab = "NDVI", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], ndvi_pix[s,], col = "red", pch=4, lwd=3)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    title(paste0("Pixel: ", pix, " (", "First Change Detected: ", dates[cps_first[pix]], ")"), cex=2)
    plot(dates, swir_pix[,1], ylab = "SWIR", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], swir_pix[s,], col = "red", pch=4, lwd=3)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, high_pix[,1], ylab = "High Albedo", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], high_pix[s,], col = "red", pch=4, lwd=3)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, low_pix[,1], ylab = "Low Albedo", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], low_pix[s,], col = "red", pch=4, lwd=3)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, soil_pix[,1], ylab = "Soil", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], soil_pix[s,], col = "red", pch=4, lwd=3)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, vege_pix[,1], ylab = "Vegetation", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], vege_pix[s,], col = "red", pch=4, lwd=3)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)

    }
    dev.off()
}


plot_all <- function(pix, data_cube, results, dates, plot_directory){
    pix <- pix
    ndvi_pix <- as.data.frame(data_cube[pix,,"ndvi"])
    swir_pix <- as.data.frame(data_cube[pix,,"swir2"])
    high_pix <- as.data.frame(data_cube[pix,,"high"])
    low_pix <- as.data.frame(data_cube[pix,,"low"])
    soil_pix <- as.data.frame(data_cube[pix,,"soil"])
    vege_pix <- as.data.frame(data_cube[pix,,"vege"])

    pdf(file = paste0(plot_directory, paste("/KR_R002-", pix, ".pdf", sep="")), width=9, height=16)
    par(mfrow=c(6,1))
    plot(dates, ndvi_pix[,1], ylab = "NDVI", xlab="Dates")
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    title(paste0("Pixel: ", pix, " (", "First Change Detected: ", dates[cps_first[pix]], ")"), cex=2)
    plot(dates, swir_pix[,1], ylab = "SWIR", xlab="Dates")
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, high_pix[,1], ylab = "High Albedo", xlab="Dates")
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, low_pix[,1], ylab = "Low Albedo", xlab="Dates")
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, soil_pix[,1], ylab = "Soil", xlab="Dates")
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, vege_pix[,1], ylab = "Vegetation", xlab="Dates")
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)

    }
    dev.off()
}

# Used for plotting rasters
get_first_cps <- function(results, dates) {
  val <- (sapply(results, function(x) {
    if (any(names(x) == "cp_info")) {
      df <- cbind(dates, 0)
      if (length(x$cp_info$cp_dates) >= 1) {
        cp_ind <- which(df[, 1] == as.numeric(min(x$cp_info$cp_dates)))[1]
        return(cp_ind)
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  }))
  return(val)
}

# Used for plotting rasters
get_cps_min_latency <- function(results, dates) {
  val <- (sapply(results, function(x) {
    if (any(names(x) == "cp_info")) {
      df <- cbind(dates, 0)
      if (length(x$cp_info$cp_dates) >= 1) {
        val1 <- min(x$cp_info$cp_detected - x$cp_info$cp_dates)
        return(val1)
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  }))
  return(val)
}

# Used for plotting rasters
get_cps_latency <- function(results, dates) {
  val <- (sapply(results, function(x) {
    if (any(names(x) == "cp_info")) {
      df <- cbind(dates, 0)
      if (length(x$cp_info$cp_dates) >= 1) {
        val1 <- (x$cp_info$cp_detected - x$cp_info$cp_dates)
        return(val1)
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  }))
  return(val)
}

# Filters for post processing
small_site_filter <- function(polys, size_thresh = 10e-5) {
  site_sizes <- expanse(polys)
  sites <- polys[site_sizes >= size_thresh]
  return(sites)
}

# Used for plotting rasters
to_polygons <- function(x, bwidth = 1e-6) {
  val <- x
  # val <- classify(val,cbind(0,NA))
  if (any(!is.na(terra::values(val)))) {
    val <- segregate(val, keep = T)
    val <- patches(val)
    val <- (fillHoles(removeDupNodes(as.polygons(val))))
    # crs(val) <- "+proj=longlat"
    # val <- buffer(val,width=bwidth)
  } else {
    val <- NULL
  }
  return(val)
}


# #############################################################################
# Plot Time Series
# #############################################################################
# KR_R002_0024
pix_13719 <- plot_all(13719, data_cube, results, dates, plots_save)
pix_14289 <- plot_all(14289, data_cube, results, dates, plots_save)
pix_14002 <- plot_all(14002, data_cube, results, dates, plots_save)
pix_17848 <- plot_all(17848, data_cube, results, dates, plots_save)
pix_23844 <- plot_all(23844, data_cube, results, dates, plots_save)

KR_R002_Pixels <- read.csv("./Pixels_to_Plot.csv")
KR_R002_Pixels <- as.data.frame(KR_R002_Pixels)
KR_R002_Pixels <- KR_R002_Pixels[,1]
for (i in KR_R002_Pixels) {
    plot_all(i, data_cube, results, dates, plots_save)
}

ts_pixels <- sample(KR_R002_Pixels, 408)
index_path <- "./pixel_index.tif"
index_rast <- rast(index_path)
samp_pts <- index_rast
keep <- c(ts_pixels)
for (ts in ts_pixels) {
    samp_pts[samp_pts == ts] <- 900
}

# #############################################################################
# Plot Rasters
# #############################################################################
library(terra)

# Load ref tif
tmp_s_path <- "./ref.tif"
tmp_s <- rast(tmp_s_path)$ref_1
values(tmp_s) <- NA

# Get first changepoints
cps_first <- get_first_cps(results, dates)
cps_latency <- get_cps_min_latency(results, dates)

# Make changepoint raster
cps_first_rast <- rast(tmp_s)
values(cps_first_rast) <- c(t(matrix(cps_first,nrow(tmp_s),ncol(tmp_s))))

dates_df <- as.data.frame(dates)
dates_df$Date_Ind <- 1:length(dates)

# Option to write raster
writeRaster(cps_first_rast, "./cps_first_priors.tif")
write.csv(keep, "./pix_sample.csv")

# #############################################################################
# Histograms of Dates
# #############################################################################
cps_shp <- readOGR("./Eval3Data/changepoints.shp")
cps_dt <- as.data.table(cps_shp)

site_24 <- cps_dt[cps_dt$Site == 24]
plot(as.Date(site_24$date), col="darkturquoise", pch=19)

site_0 <- cps_dt[cps_dt$Site == 0]
plot(as.Date(site_0$date), col="darkturquoise", pch=19)

site_2 <- cps_dt[cps_dt$Site == 2]
plot(as.Date(site_2$date), col="darkturquoise", pch=19)

site_3 <- cps_dt[cps_dt$Site == 3]
plot(as.Date(site_3$date), col="darkturquoise", pch=19)

site_5 <- cps_dt[cps_dt$Site == 5]
plot(as.Date(site_5$date), col="darkturquoise", pch=19)


for (i in 1:length(dates)) {
    print(paste0("From ", dates[i], " to ", dates[i+1], ": ", 
    difftime(dates[i], dates[i+1], units="days"), " days"))
}

