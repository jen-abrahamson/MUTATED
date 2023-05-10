#*******************************************************************************
# Project: MUTATED - roboBayes
# Script Purpose: test temporal window for pseudo online mode
# Date: 2022-08-01
# Author: Jenna Abrahamson
#*******************************************************************************
# Load libraries
library(terra)
library(fields)
library(inlmisc)
library(optparse)
library(RJSONIO)
library(data.table)
library(geojsonio)
library(rgdal)

# Source heuristic filter
source("./bas/robobayes/functions/post_filters.R")

# load results data
load("./MUTATED/data/BH_R001_ACC/processed/BH_R001_ACC_Tune8_L8.Rdata")

# #############################################################################
# Plotting functions
# #############################################################################

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
# Functions to construct data paths by region
# #############################################################################


define_paths <- function(region, local) {
    tmp_s_path <- paste0("/rsstu/users/j/jmgray2/SEAL/MUTATED/data/", region, "_ACC/interim/ref.tif")
    annotations_path <- paste0("/rsstu/users/j/jmgray2/SEAL/MUTATED/data/annotations/region_models/", 
                               region, ".geojson")
    proj4_path <- paste0("/rsstu/users/j/jmgray2/SEAL/MUTATED/data/", region, "_ACC/interim/proj4.txt")
    if (local == TRUE) {
        tmp_s_path <- paste0("/Volumes", tmp_s_path)
        annotations_path <- paste0("/Volumes", annotations_path)
        proj4_path <- paste0("/Volumes", proj4_path)
    }
    files_list <- list("tmp_s" = tmp_s_path, "annotations" = annotations_path, "proj4" = proj4_path)
    return(files_list)
}


paths <- define_paths("BH_R001", local=TRUE)




# #############################################################################
# Plot Rasters (optional)
# #############################################################################

# Load ref tif
tmp_s <- rast(paths$tmp_s)$ref_1
values(tmp_s) <- NA

# Get first changepoints
cps_first <- get_first_cps(results, dates)
cps_latency <- get_cps_min_latency(results, dates)

# Make changepoint raster
cps_first_rast <- rast(tmp_s)
values(cps_first_rast) <- c(t(matrix(cps_first,nrow(tmp_s),ncol(tmp_s))))


# #############################################################################
# Running Post-Processs
# #############################################################################
option_list <- list(
  make_option(c("--params"),
    type = "character",
    default = "../../../cirrus/robobayes.json"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Need example params file, change if necessary
params <- fromJSON("./BH_R001_T8.json")$ROBOBAYES
out_tag <- params$OUT_TAG
region <- params$REGION
data_dir <- params$DATA_DIR
bands <- params$BANDS 
filt_1 <- params$FILT_NDVI_AMP
filt_2 <- params$FILT_NDVI_MEAN
filt_3 <- params$FILT_VEGE_AMP
filt_4 <- params$FILT_LA_VAR
filt_5 <- params$FILT_NDVI_CURRAMP_LIM
filt_6 <- params$FILT_NDVI_AMP_LIM
filt_7 <- params$FILT_HA_LIM



######## Set up spatial information

# Load ref tif
tmp_s_path <- paths$tmp_s
tmp_s <- rast(tmp_s_path)$ref_1
values(tmp_s) <- NA

# Load annotations
annotations_path <- paths$annotations
annotations <- geojson_read(annotations_path, what = "sp") 
# Load coordinate system for ref tif
proj4_path <- paths$proj4
ref_crs <- read.table(proj4_path, header=F, sep="|", stringsAsFactors = F)[1, "V1"] 
annotations <- spTransform(annotations, ref_crs) 

# Load base
base_tif <- rast("/Volumes/rsstu/users/j/jmgray2/SEAL/JennaAbrahamson/MUTATED/pseudo_online/gif_plots/BH_R001/base.tif")
ext(base_tif) <- ext(annotations)

# Function to create pseudo online plots for animations
run_pseudo_online <- function(window, dates, results, tmp_s, annotations, base_tif){
    # Get cps_first raster and dates list
    cps_first <- get_first_cps(results, dates)

    # Subset annotations for plotting
    pos_annotations <- subset(annotations, status=="positive_annotated")
    other_annotations <- subset(annotations, status!="positive_annotated")

    # Create empty spat vector to append to
    polys_final <- vect()

    # Set up moving window params
    date_list <- c(1:length(dates))
    w_floor <- floor(window / 2) # number of obs to consider before/after
    n_obs <- max(date_list)

    # Run moving temporal window
    for (i in (w_floor + 1):(length(dates) - w_floor)) {

        # Determine observation of interest
        window_obs <- date_list[(i - w_floor):(i + w_floor)]
        cps <- cps_first
        cps[cps < min(window_obs) | cps > max(window_obs)] <- NA

        # Run filtering
        cps_first_filtered <- param_filter(results, cps, filt_ndvi_amp=0.15, 
                                        filt_ndvi_mean=0.4, filt_vege_amp=0.15, 
                                        filt_la_var=0.05, filt_ndvi_curramp_lim=0.02, 
                                        filt_ndvi_amp_lim=0.02, filt_ha_lim=0.4)

        # Get results into raster
        values(tmp_s) <- NA
        cps_first_rast_filtered <- tmp_s
        values(cps_first_rast_filtered) <- c(t(matrix(cps_first_filtered,nrow(tmp_s),ncol(tmp_s))))
    
        # Print as a test check
        print(summary(cps))

        # Load annotations adn subset
        pos_annot <- subset(pos_annotations, start_date <= dates[i])
        other_annot <- subset(other_annotations, start_date <= dates[i])
        # Get into polygons
        polys <- to_polygons(cps_first_rast_filtered)

        if (length(polys) == 0) {

            next
        } else {
        polys_filtered <- (small_site_filter(polys, size_thresh = 9000))
        }
        if (length(polys_filtered) == 0) {
            next
        } else {
            if (nchar(i) == 1) {
                file_name = paste0(plot_dir, "000", i, ".png", sep="")
            }
            if (nchar(i) == 2) {
                file_name = paste0(plot_dir, "00", i, ".png", sep="")
            }
            if (nchar(i) == 3) {
                file_name = paste0(plot_dir, "0", i, ".png", sep="")
            }
            png(file = file_name, 
                width=875, height=950, units="px", bg="white", pointsize=12)
            plotRGB(base, alpha=170)
            polys_final <- rbind(polys_final, polys_filtered)
            plot(polys_final, col="deepskyblue4", add=T, border="black")
            plot(polys_filtered, col = "cyan2", border = "black", add=T)
            plot(pos_annot, add=T, border="red", col=NA, lwd=2.75)
            plot(other_annot, add=T, border="darkorchid3", col=NA, lwd=2.75)
            text(459999, 2901419, paste0(dates[min(window_obs)], " to ", dates[max(window_obs)], sep=""), 
                col="black", cex=1.85, lwd=2, font=2)
            dev.off()
        }

    }
}

# Set plotting dir and run function
plot_dir <- "./gif_plots/BH_R001/window_10/"
test <- run_pseudo_online(window=10, dates, results, tmp_s, annotations, base)


# Command line code for creating gifs (Don't run here)
convert -delay 20 -loop 0 *.png KR_R002_window10.gif 
mogrify -layers 'optimize' -fuzz 7% US_R004_Signals.gif



