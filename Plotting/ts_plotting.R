#*******************************************************************************
# Project: Plot Time Series
# Script Purpose: Code and functions for easily plotting time series
# Date: 2022-02-22
# Author: Jenna Abrahamson
#*******************************************************************************

#***********************************
# Set user inputs/load data
#***********************************
# Change these to desired region data

# Load original Rdata
og_data <- load("/Volumes/rsstu/users/j/jmgray2/SEAL/MUTATED/data/BH_R001_ACC/BH_R001_ACC.Rdata")
# Load results data
results_data <- load("/Volumes/rsstu/users/j/jmgray2/SEAL/MUTATED/data/BH_R001_ACC/processed/BH_R001_ACC_BH_SpikeReturn.Rdata")

#***********************************
# Clean/order all data
#***********************************
# Run everything from here until next subsection

sort_ind <- sort.int(hls_dates, index.return = T)$ix
hls_sorted <- sort(hls_dates)
all_t <- as.integer(hls_sorted)

scale_factor <- 1e-4

# reorder data to match reordered time
qa_df_sorted <- as.matrix(qa_df)[, sort_ind]
ndvi_df_sorted <- as.matrix(ndvi_df)[, sort_ind]
swir2_df_sorted <- as.matrix(swir2_df)[, sort_ind] * scale_factor
high_df_sorted <- as.matrix(high_df)[, sort_ind]
low_df_sorted <- as.matrix(low_df)[, sort_ind]
soil_df_sorted <- as.matrix(soil_df)[, sort_ind]
veg_df_sorted <- as.matrix(veg_df)[, sort_ind]

# identify additional bad points
remove_ind <- (swir2_df_sorted <= 0) | (ndvi_df_sorted < 0) |
    (ndvi_df_sorted > 1) | is.nan(swir2_df_sorted) |
    is.nan(ndvi_df_sorted)

remove_ind <- remove_ind | (qa_df_sorted != 0)

# remove them
ndvi_df_sorted[remove_ind] <- NA
swir2_df_sorted[remove_ind] <- NA
high_df_sorted[remove_ind] <- NA
low_df_sorted[remove_ind] <- NA
soil_df_sorted[remove_ind] <- NA
veg_df_sorted[remove_ind] <- NA

normalize <- function(x) (x - 0) / (10000 - 0)

#***********************************
# Plotting NDVI & SWIR time series with detected change dates
#***********************************

# Plotting function
plot_NDVI_SWIR <- function (pixel, ndvi_df, swir_df, results_list) {
    ndvi_pixel <- ndvi_df[pixel,]
    swir_pixel <- swir_df[pixel,]
    cp_dates <- results_list[pixel][[1]]$cp_info$cp_dates
    par(mfrow=c(2,1))
    plot(hls_sorted, ndvi_df[pixel,], xlab = "Date", ylab="NDVI", ylim=c(0,1))
    for (d in cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    legend("topright", legend="Detected Change Date", lty= 1, lwd = 2, col = "red")
    plot(hls_sorted, swir_df[pixel,], xlab = "Date", ylab="SWIR", ylim=c(0,1))
    for (d in cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    legend("topright", legend="Detected Change Date", lty= 1, lwd = 2, col = "red")
}

# +++++++++++++++++
# Example 
plot_NDVI_SWIR(79, ndvi_df_sorted, swir2_df_sorted, results)

#***********************************
# Plotting NDVI & SWIR time series with detected change dates and outliers
#***********************************
# Shows plots for NDVI and SWIR with all dtected change dates, and outlier
# points removed by the spike filter

# Plotting function
plot_NDVI_SWIR <- function (pixel, ndvi_df, swir_df, results_list) {
    ndvi_pixel <- ndvi_df[pixel,]
    swir_pixel <- swir_df[pixel,]
    cp_dates <- results_list[pixel][[1]]$cp_info$cp_dates
    pixel_spikes <- c(spikes[pixel,])
    unlist_spikes <- unlist(pixel_spikes)
    par(mfrow=c(2,1))
    plot(hls_sorted, ndvi_df[pixel,], xlab = "Date", ylab="NDVI", ylim=c(0,0.8))
    for (d in cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    for (s in unlist_spikes){
        points(hls_sorted[s], ndvi_df[pixel,s], col = "blue", pch=4, lwd = 3, cex=1.1)
    }
    legend("topright", legend=c("Removed Outlier", "Detected Change Date"), pch=c(4, NA), lty= c(NA, 1), lwd = c(3,2), col = c("blue", "red"), bty=)

    plot(hls_sorted, swir_df[pixel,], xlab = "Date", ylab="SWIR", ylim=c(0,0.8))
    for (d in cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    for (s in unlist_spikes){
    points(hls_sorted[s], swir_df[pixel,s], col = "blue", pch=4, lwd = 3, cex=1.1)
    }
    legend("topright", legend=c("Removed Outlier", "Detected Change Date"), pch=c(4, NA), lty= c(NA,1), lwd = c(3,2), col = c("blue", "red"))
}

# +++++++++++++++++
# Example 
plot_NDVI_SWIR(3648, ndvi_df_sorted, swir2_df_sorted, results)

cps <- get_first_cps(results, dates)




# Testing with data_cube
plot_NDVI_SWIR(33503)


pix <- 35161
ndvi_pix <- as.data.frame(data_cube[pix,,"ndvi"])
swir_pix <- as.data.frame(data_cube[pix,,"swir2"])
high_pix <- as.data.frame(data_cube[pix,,"high"])
low_pix <- as.data.frame(data_cube[pix,,"low"])
soil_pix <- as.data.frame(data_cube[pix,,"soil"])
vege_pix <- as.data.frame(data_cube[pix,,"vege"])

par(mfrow=c(6,1))
plot(dates, ndvi_pix[,1])
for (s in spikes[pix]){
    points(dates[s], ndvi_pix[s,], col = "red", pch=4, lwd=4)
}
for (d in results[[pix]]$cp_info$cp_dates) {
    print(d) # Check dates
    abline(v = d, col = "red", lwd=2)
}
title(pix)
plot(dates, swir_pix[,1])
for (s in spikes[pix]){
    points(dates[s], swir_pix[s,], col = "red", pch=4, lwd=4)
}
for (d in results[[pix]]$cp_info$cp_dates) {
    print(d) # Check dates
    abline(v = d, col = "red", lwd=2)
}
plot(dates, high_pix[,1])
for (s in spikes[pix]){
    points(dates[s], high_pix[s,], col = "red", pch=4, lwd=4)
}
for (d in results[[pix]]$cp_info$cp_dates) {
    print(d) # Check dates
    abline(v = d, col = "red", lwd=2)
}
plot(dates, low_pix[,1])
for (s in spikes[pix]){
    points(dates[s], low_pix[s,], col = "red", pch=4, lwd=4)
}
for (d in results[[pix]]$cp_info$cp_dates) {
    print(d) # Check dates
    abline(v = d, col = "red", lwd=2)
}
plot(dates, soil_pix[,1])
for (s in spikes[pix]){
    points(dates[s], soil_pix[s,], col = "red", pch=4, lwd=4)
}
for (d in results[[pix]]$cp_info$cp_dates) {
    print(d) # Check dates
    abline(v = d, col = "red", lwd=2)
}
plot(dates, vege_pix[,1])
for (s in spikes[pix]){
    points(dates[s], vege_pix[s,], col = "red", pch=4, lwd=4)
}
for (d in results[[pix]]$cp_info$cp_dates) {
    print(d) # Check dates
    abline(v = d, col = "red", lwd=2)
}


plotAll <- function(pix, data_cube, results, dates){
    pix <- pix
    ndvi_pix <- as.data.frame(data_cube[pix,,"ndvi"])
    swir_pix <- as.data.frame(data_cube[pix,,"swir2"])
    high_pix <- as.data.frame(data_cube[pix,,"high"])
    low_pix <- as.data.frame(data_cube[pix,,"low"])
    soil_pix <- as.data.frame(data_cube[pix,,"soil"])
    vege_pix <- as.data.frame(data_cube[pix,,"vege"])

    par(mfrow=c(6,1))
    plot(dates, ndvi_pix[,1], ylab = "NDVI", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], ndvi_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    title(paste0("Pixel: ", pix))
    plot(dates, swir_pix[,1], ylab = "SWIR", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], swir_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, high_pix[,1], ylab = "High Albedo", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], high_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, low_pix[,1], ylab = "Low Albedo", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], low_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, soil_pix[,1], ylab = "Soil", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], soil_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    plot(dates, vege_pix[,1], ylab = "Vegetation", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], vege_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
}


plotAll(533, data_cube, results, dates)

png(file = paste0(plot_directory, paste(cell, ylab, ".png", sep="")), width=1200, height=400)
      plot(hls_sorted, pix, main = paste0("Pixel:", cell), xlab = "Date", ylab=ylab, ylim= c(0,1))
      points(hls_sorted[spikes], pix[spikes], col = "red", pch = 4,lwd = 4, cex=1.5)
      dev.off()


saveAll <- function(pix, data_cube, results, dates, plot_directory){
    pix <- pix
    ndvi_pix <- as.data.frame(data_cube[pix,,"ndvi"])
    swir_pix <- as.data.frame(data_cube[pix,,"swir2"])
    high_pix <- as.data.frame(data_cube[pix,,"high"])
    low_pix <- as.data.frame(data_cube[pix,,"low"])
    soil_pix <- as.data.frame(data_cube[pix,,"soil"])
    vege_pix <- as.data.frame(data_cube[pix,,"vege"])

    png(file = paste0(plot_directory, paste("KR_R002-", pix, ".png", sep="")), width=500, height=250, units="mm", res=500)
    par(mfrow=c(3,2))
    plot(dates, ndvi_pix[,1], ylab = "NDVI", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], ndvi_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    title(paste0("Pixel: ", pix))
    plot(dates, swir_pix[,1], ylab = "SWIR", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], swir_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    title(paste0("Pixel: ", pix))
    plot(dates, high_pix[,1], ylab = "High Albedo", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], high_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    title(paste0("Pixel: ", pix))
    plot(dates, low_pix[,1], ylab = "Low Albedo", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], low_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    title(paste0("Pixel: ", pix))
    plot(dates, soil_pix[,1], ylab = "Soil", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], soil_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    title(paste0("Pixel: ", pix))
    plot(dates, vege_pix[,1], ylab = "Vegetation", xlab="Dates")
    for (s in spikes[pix]){
        points(dates[s], vege_pix[s,], col = "red", pch=4, lwd=4)
    }
    for (d in results[[pix]]$cp_info$cp_dates) {
        print(d) # Check dates
        abline(v = d, col = "red", lwd=2)
    }
    title(paste0("Pixel: ", pix))
    dev.off()
}

plot_directory <- "/Volumes/rsstu/users/j/jmgray2/SEAL/JennaAbrahamson/MUTATED/Outputs/Updated_roboBayes/KR_R002_Plots/"

saveAll(533, data_cube, results, dates, plot_directory)

# Save a bunch from Korea
library(sf)
library(rgdal)
test_pts <- st_read("/Volumes/rsstu/users/j/jmgray2/SEAL/JennaAbrahamson/MUTATED/Outputs/OutlierMethodsTest/KR_PixelsCheck.shp")
pix_to_plot <- test_pts$Pixel_No

for (pix in pix_to_plot) {
    saveAll(pix, data_cube, results, dates, plot_directory)
}

for (pix in pix_to_plot) {
    print(results[[pix]]$cp_info$cp_dates)
}