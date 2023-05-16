#*******************************************************************************
# Project: MUTATED - roboBayes
# Script Purpose: Functions for building tree features
# Date Created: 2022-03-16
# Author: Jenna Abrahamson
#*******************************************************************************
# Load libraries
library(terra)
library(data.table)
library(matrixStats)
library(glcm)

#' relative_change - calculate relative change between values
#' 
#'@description
#' Takes in an original value and a new value and calculates 
#' the relative change.
#' 
#' @param final_val The new (or final) value (numeric)
#' @param og_val The original (or initial) value (numeric)
#' 
#' @return relative change (numeric)
#' 
#' @export
relative_change <- function(og_val, final_val) {
    x <- (final_val - og_val) / og_val
    return(x)
}


#' summarize_signal - summarize values from a specific signal
#' 
#'@description
#' Takes in the signal values from the data cube and generates the 
#' mediain, variance, and scene median.
#' 
#' @param signal the name of the signal (character)
#' @param data_cube data cube object (array)
#' 
#' @return summarized signal data (data.frame)
#' 
#' @export
# helper function to summarize signal
summarize_signal <- function(signal, data_cube) {
    # Isolate single signal from data cube
    sig <- as.data.table(as.data.frame(data_cube[,,signal]))

    # Get list of dates
    date_select <- colnames(sig)

    # Calculate summary stats across all obs
    sig <- sig[, paste0(signal, "_median") := rowMedians(as.matrix(.SD), na.rm=T), .SDcols = date_select]
    sig <- sig[, paste0(signal, "_var") := rowVars(as.matrix(.SD), na.rm=T), .SDcols = date_select]
    sig <- sig[, paste0(signal, "_mean") := rowMeans(as.matrix(.SD), na.rm=T), .SDcols = date_select]
    sig <- sig[, paste0(signal, "_scene_median") := median(as.matrix(sig[,1:(dim(sig)[2]-2)]), na.rm=T)]    

    # Turn 0s to -> 0.0000001 to avoid inf values
    sig[sig == 0] <- 0.0000001
    # Calculate relative change of pixel val compared to scene median
    df <- as.data.frame(sig)
    rel_median <- relative_change(df[,ncol(df)], df[,(ncol(df)-3)])
    sig <- cbind(sig, rel_median)
    names(sig)[names(sig) == "rel_median"] <- paste0(signal, "_rel_change")
    sig$index <- 1:dim(sig)[1]
    return(as.data.frame(sig))
}

#' calc_texture_metrics - calculate texture features
#' 
#'@description
#' Takes in the signal train data.frame and generates
#' texture features
#' 
#' @param sig_train_dt signal data returned from summarize_signal (data.frame)
#' @param tmp_s reference tif (SpatRaster)
#' @param rast_index A raster containing the pixel indices (SpatRaster)
#' 
#' @return texture metrics for specified signal (data.frame)
#' 
#' @export
calc_texture_metrics <- function(sig_feat, ref_tif, rast_index) {
    # Initialize df of pixel indices
    tex_feat <- data.frame(matrix(ncol=1, nrow=(dim(sig_feat)[1])))
    tex_feat$Pixel_Index <- values(rast_index)
    tex_feat <- tex_feat[,2]
    colnames(tex_feat) <- "Pixel_Index"

    # Loop through all dates
    for (date in colnames(sig_feat)[1:(length(colnames(sig_feat))-6)]){
        sig_feat <- as.data.frame(sig_feat)
        sig_data <- as.data.table(sig_feat[, date])
        sig_data$V1[is.infinite(sig_data$V1)] <- NA # Ensure there's no infs

        # Catch dates with mostly NA values and skip
        if (sum(is.na(sig_data)) > (dim(sig_feat)[1] * 0.95)){
            print(paste0(date, " mostly NA values, skipping!"))
            next
        }
        if (var(sig_data, na.rm=TRUE) == 0) {
            print(paste0(date, " no variance to signal values, skipping!"))
            next
        } else {
            # Initialize signal raster
            colnames(sig_data) <- "sig_data"
            sig_rast <- ref_tif
            values(sig_rast) <- NA
            values(sig_rast) <- c(t(matrix(c(sig_data$sig_data),nrow(ref_tif),ncol(ref_tif))))

            # Calculate texture metrics
            library(raster) # 
            texture_metrics <- glcm(raster(sig_rast), 
                                window = c(3,3), # use a 3 x 3 window
                                shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)), 
                                na_opt = "ignore", # ignore any NA values
                                statistics = c("mean",
                                               "variance",
                                               "homogeneity",
                                               "contrast",
                                               "entropy", 
                                               "dissimilarity",
                                               "second_moment", 
                                               "correlation"))

            # Find edges
            hpf <- matrix(c(-1, -1, -1, -1, 9, -1, -1, -1, -1), nrow=3)
            hpf_filter <- focal(raster(sig_rast), w = hpf)

            # Create colnames to use
            edge_col <- paste0("hpf_", date)
            mean_col <- paste0("mean_", date)
            var_col <- paste0("var_", date)
            homo_col <- paste0("homo_", date)
            cont_col <- paste0("cont_", date)
            entr_col <- paste0("entr_", date)
            diss_col <- paste0("diss_", date)
            sec_col <- paste0("sec_", date)
            corr_col <- paste0("corr_", date)
            tex_names <- c(edge_col, mean_col, var_col, homo_col, cont_col, entr_col, diss_col, sec_col, corr_col)
        
            # Add to sig_data data.table
            sig_data$edges <- values(hpf_filter)
            sig_data$glcm_mean <- values(texture_metrics$glcm_mean)
            sig_data$glcm_var <- values(texture_metrics$glcm_variance)
            sig_data$glcm_homo <- values(texture_metrics$glcm_homogeneity)
            sig_data$glcm_cont <- values(texture_metrics$glcm_contrast)
            sig_data$glcm_entr <- values(texture_metrics$glcm_entropy)
            sig_data$glcm_diss <- values(texture_metrics$glcm_dissimilarity)
            sig_data$glcm_sec <- values(texture_metrics$glcm_second_moment)
            sig_data$glcm_corr <- values(texture_metrics$glcm_correlation)

            # colnames
            colnames(sig_data)[(dim(sig_data)[2] - 8):dim(sig_data)[2]] <- tex_names
            tex_feat <- cbind(tex_feat, sig_data)
            print(paste0("Finished Extracting ", date, "!"))
        }
    }
    dates_list <- colnames(tex_feat)
    tex_feat <- as.data.frame(tex_feat)
    
    hpf_text <- as.data.table(tex_feat[, c(dates_list[grep('hpf', dates_list)])])
    hpf_text <- hpf_text[, hpf_text := rowMedians(as.matrix(.SD), na.rm=T)]
    mean_text <- as.data.table(tex_feat[, c(dates_list[grep('mean', dates_list)])])
    mean_text <- mean_text[, mean_text := rowMedians(as.matrix(.SD), na.rm=T)]
    var_text <- as.data.table(tex_feat[, c(dates_list[grep('var', dates_list)])])
    var_text <- var_text[, var_text := rowMedians(as.matrix(.SD), na.rm=T)]
    homo_text <- as.data.table(tex_feat[, c(dates_list[grep('homo', dates_list)])])
    homo_text <- homo_text[, homo_text := rowMedians(as.matrix(.SD), na.rm=T)]
    cont_text <- as.data.table(tex_feat[, c(dates_list[grep('cont', dates_list)])])
    cont_text <- cont_text[, cont_text := rowMedians(as.matrix(.SD), na.rm=T)]
    entr_text <- as.data.table(tex_feat[, c(dates_list[grep('entr', dates_list)])])
    entr_text <- entr_text[, entr_text := rowMedians(as.matrix(.SD), na.rm=T)]
    diss_text <- as.data.table(tex_feat[, c(dates_list[grep('diss', dates_list)])])
    diss_text <- diss_text[, diss_text := rowMedians(as.matrix(.SD), na.rm=T)]
    sec_text <- as.data.table(tex_feat[, c(dates_list[grep('sec', dates_list)])])
    sec_text <- sec_text[, sec_text := rowMedians(as.matrix(.SD), na.rm=T)]
    corr_text <- as.data.table(tex_feat[, c(dates_list[grep('corr', dates_list)])])
    corr_text <- corr_text[, corr_text := rowMedians(as.matrix(.SD), na.rm=T)]

    # Combine to one data.frame
    final_texture <- as.data.frame(cbind(tex_feat$Pixel_Index, hpf_text$hpf_text, 
                                         mean_text$mean_text, var_text$var_text, homo_text$homo_text, 
                                         cont_text$cont_text, entr_text$entr_text, diss_text$diss_text, 
                                         sec_text$sec_text, corr_text$corr_text))
    
    colnames(final_texture) <- c("pixel_index", "glcm_hpf", "glcm_mean", "glcm_var", "glcm_homo", "glcm_cont", "glcm_entr", "glcm_diss", "glcm_sec", "glcm_corr")
    final_texture <- data.frame(sapply(final_texture, function(x) ifelse(is.nan(x), NA, x)))
    return(final_texture)
}


#' build_input_features - generate input features
#' 
#'@description
#' Takes in data and change polygons and generates input features
#' to run random forest prediciton
#' 
#' @param data_cube data cube object (array)
#' @param change_polys detected change sites (geojson, SpatVector)
#' @param ref_tif reference tif (SpatRaster)
#' @param possible_signals vector of possible signals from params json (vector)
#' 
#' @return all input data (data.frame)
#' 
#' @export
# Build input features
build_input_features <- function(data_cube, results, ref_tif, possible_signals) {
    # Pre-process data cube
    data_cube <- input_data_check(data_cube, dates, index)
    # scale data to be between 0 and 1
    data_cube <- signal_scaling(data_cube, possible_signals)
    # Elimanate data that does not obey the data boundsß
    data_cube <- domain_filter(data_cube, possible_signals)
    # Eliminate data contaminated by clouds, aerosols, etc.
    data_cube <- quality_filter(data_cube)

    # Calculate brightness and add to data.cube
    data_cube[,,"brightness"] <- sqrt((data_cube[,,"red"] * data_cube[,,"red"])/
                             (data_cube[,,"green"] * data_cube[,,"green"])/2)

    # Take last 10 observations
    data_cube <- data_cube[,(dim(data_cube)[2] - 9):dim(data_cube)[2],]

    # Generate cps_first raster from results
    cps_first <- get_first_cps(results, dates)
    cps_first_rast <- ref_tif
    values(cps_first_rast) <- NA
    values(cps_first_rast) <- c(t(matrix(cps_first,nrow(ref_tif),ncol(ref_tif))))

    # Create raster of pixel indices
    rast_index <- ref_tif
    values(rast_index) <- NA
    values(rast_index) <- c(t(matrix(c(1:length(values(ref_tif))),nrow(ref_tif),ncol(ref_tif))))
    print('get clean data.frame of change pixel features')
    # Get clean data.frame of change pixel indices
    cps_first_rast[is.nan(cps_first_rast)] <- NA # Ensure there's no NaNs
    change_pixels <- c(rast_index, cps_first_rast)
    names(change_pixels) <- c("index", "cp")

    # Only keep indices that have an actual cp
    change_pixels$index[is.na(change_pixels$cp)] <- NA
    change_df <- as.data.frame(change_pixels$index)

    # ~~~~~~~~~~~~~~~~~~~~~~
    # Extract spectral data
    # ~~~~~~~~~~~~~~~~~~~~~~
    # Initiate data.table to add feature data to
    features <- as.data.table(change_df)
    colnames(features) <- c("Pixel_Index")
    features <- na.omit(features)

    # Summarize all spectral signals and append to features data.table
    signals <- c("red", "green", "blue", "nir", "swir1", "swir2", "ndvi", "brightness")
    for (sig in signals) {
        sig_feat <- summarize_signal(sig, data_cube)
        features <- cbind(features, sig_feat[features$Pixel_Index, c(11:15)])

        # Calc texture metrics
        sig_text <- calc_texture_metrics(sig_feat, tmp_s, rast_index)
        sig_text <- sig_text[order(sig_text$pixel_index),]

        # Create labels
        texture_labels <- colnames(sig_text)[2:dim(sig_text)[2]]
        texture_labels <- paste0(sig, "_", texture_labels)
        colnames(sig_text)[2:dim(sig_text)[2]] <- texture_labels
        features <- cbind(features, sig_text[features$Pixel_Index, c(2:dim(sig_text)[2])])
    }

    # Put into final data.table
    features <- as.data.table(features)

    # Replace all inf values with NA
    is.na(features) <- do.call(cbind, lapply(features, is.infinite))

    # Impute NA valuesß
    for (col in 2:ncol(features)){
        features <- as.data.frame(features)
        features[,col][is.na(features[,col])] <- median(features[,col],na.rm=TRUE)
    }

    # return the final file with all training data
    return(features)
}

#' filter_results - filter detect change pixels
#' 
#'@description
#' Takes in RF results and filters detected change sites
#' 
#' @param features data.frame of inputs (data.frame)
#' @param rf_model random forest model (randomForest obj)
#' @param cps_first detected change sites, with start date (vector)
#' @param cps_last detected change sites, with end date (vector)
#' @param cps_cps detected change sites, binary change (vector)
#' @param ref_tif reference tif (SpatRaster)
#' @param resutls roboBayes results object (list)
#' 
#' @return construction change sites (SpatVector)
#' 
#' @export

filter_results <- function(features, rf_model, cps_first, cps_last, cps, ref_tif) {
    # Filter any inf values
    is.na(features) <- do.call(cbind, lapply(features, is.infinite))
    features <- as.data.frame(features)

    # Predict with RF model
    pred <- predict(rf_model, features[,2:dim(features)[2]])

    # Create data table to screen with
    filter <- as.data.table(features$Pixel_Index)
    filter$Prediction <- as.character(pred)
    false_pos <- filter$V1[filter$Prediction == "Other"]

    # Screen vectors
    cps_first[false_pos] <- NA
    cps_last[false_pos] <- NA 
    cps[false_pos] <- 0

    # Convert to rasters
    cps_first_rast <- ref_tif
    values(cps_first_rast) <- c(t(matrix(cps_first,nrow(ref_tif),ncol(ref_tif))))

    cps_last_rast <- ref_tif
    values(cps_last_rast) <- c(t(matrix(cps_last,nrow(ref_tif),ncol(ref_tif))))

    cps_rast <- ref_tif
    values(cps_rast) <- c(t(matrix(cps,nrow(ref_tif),ncol(ref_tif))))

    # Combine into one SpatVector
    tifs <- c(cps_rast, cps_first_rast, cps_last_rast)
    print(summary(tifs))
    names(tifs) <- c("cps_rast", "cps_first_rast", "cps_last_rast")

    return(tifs)
}
