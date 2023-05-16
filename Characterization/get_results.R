#*******************************************************************************
# Project: MUTATED - roboBayes
# Script Purpose: Generate Results w/ Tree Filter
# Date: 2023-03-04
# Author: Jenna Abrahamson
#*******************************************************************************
Sys.setenv(PROJ_LIB="/usr/local/usrapps/jmgray2/jnabraha/bocpd/bin/proj")
library(terra)
library(fields)
library(inlmisc)
library(optparse)
library(RJSONIO)
library(data.table)
library(xgboost)
source("./data_prep.R")
source("./build_features.R")
source("./post_extract_quantities.R")

# Load possible signals
possible_signals <- c("blue", "green", "red", "nir", "swir1", "swir2",
    "ndvi", "high", "low", "vege", "soil", "qa", "brightness")

# Set data dir
data_dir <- 
output_dir <- 
model_src <- "./globalXG.Rdata"

# Create list of file paths to use
regions <- list.files(data_dir)
regions <- regions[grep("R00", regions)]

# Loop through regions and generate results
for (region in regions) {
    # Reload terra
    library(terra)

    # Load og data cube
    og_src <- paste0(data_dir, region, "/")
    sub_dirs <- list.files(og_src)
    sub_dirs <- sub_dirs[grepl("^[[:digit:]]+", sub_dirs)]

    # Loop through mgrs tiles
    for (sub in sub_dirs) {
        # Load og data cube
        sub_src <- paste0(data_dir, region, "/", sub, "/bas_intermediates/ACC.Rdata")
        if (file.exists(sub_src) == FALSE) {
            next
        }
        load(sub_src)

        # Load data cube results
        sub_results <- paste0(data_dir, region, "/", sub, "/bas_intermediates/interim/ACC.Rdata")
        if (file.exists(sub_results) == FALSE) {
            next
        }
        load(sub_results)
    
        # Load in ref tif
        tmp_s_src <- paste0(data_dir, region, "/", sub, "/bas_intermediates/interim/tmp.tif")
        if (file.exists(tmp_s_src) == FALSE) {
            next
        }

        # Go through post-processing
        tmp_s <- rast(tmp_s_src)
        tmp_s <- tmp_s[names(tmp_s)[1]]
        terra::values(tmp_s) <- NA

        unfiltered_cps_first <- get_first_cps(results, dates)

        features <- build_input_features(data_cube, 
                                 results, 
                                 tmp_s, 
                                 possible_signals)


        # Load in RF model
        print("Made it to loading RF model")
        model <- get(load(model_src))

        # Run RF filter
        filtered_polys <- filter_results(features, 
                                        model, 
                                        unfiltered_cps_first, 
                                        tmp_s)

        n_polys <- as.numeric(nrow(filtered_polys))
        if (isTRUE(n_polys > 0) == FALSE) {
            print("No polygons found!")
            next
        }
        
        out_file <- paste0(output_dir, region, "_", sub, ".geojson")

        # Output the final geojson
        print("check_final")
        
        # Write out polys
        writeVector(filtered_polys, file = out_file, filetype = "geojson", overwrite=T)
        print(paste0("Finished running ", region, "_", sub, "!"))
    }
}

