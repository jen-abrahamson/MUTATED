#*******************************************************************************
# Project: roboBayes - Spike Filter Grid Search
# Script Purpose: Find optimal spike filter parameters for each signal 
#                 (central moving window spike filter)
# Date: 2022-03-06
# Author: Jenna Abrahamson
#*******************************************************************************

# Load necessary libraries
library(paramtest)
library(data.table)

#***********************************
# Set file paths
#***********************************
# Load original Rdata file
original_Rdata <- "./data/KR_R002/KR_R002.Rdata"
load(original_Rdata) 

# Set path for to plot directory if desired
plot_directory <- "./Spike_GridSearch/TS_Plots/"

# Set paths to training data
ndvi_train <- "./Spike_GridSearch/NDVI_Train.csv"
swir_train <- "./Spike_GridSearch/SWIR2_Train.csv"
ha_train <- "./Spike_GridSearch/HA_Train.csv"
la_train <- "./Spike_GridSearch/LA_Train.csv"
soil_train <- "./Spike_GridSearch/Soil_Train.csv"
veg_train <- "./Spike_GridSearch/Veg_Train.csv"

#*******************************************************************************
# Functions
#*******************************************************************************

#***********************************
# Main spike filter function - spike_filter()
#***********************************
#' Detect outlier spikes in a time series
#' 
#'@description
#' Takes in a single pixel's time series and uses a moving window to flag
#' spikes that are outside of a set standard deviation threshold, or have
#' a difference in spike amplitude that violates a set threshold.
#' 
#' @param signal The time series from which to flag spikes (data.frame)
#' @param window The length of the moving window (integer)
#' @param threshold The # of std away from median to flag spike (numeric)
#' @param spike_amp Spike amplitude (units of value being compared) (numeric)
#' @param timeframe Discard if observations span more than x days (integer)
#' 
#' @return A vector of pixel indices to be screened from the input time series
#' @export

spike_filter <- function(signal, window=11, threshold=0.1, spike_amp=0.01, timeframe=100, dates=hls_sorted) {
  sig_df <- as.data.table(signal)
  pixel <- cbind(sig_df, "timestamp" = 1:length(dates))
  pixel <- na.omit(pixel)
  w_floor <- floor(window / 2)
  pixel$spike <- 0
  spike_dates <- c()
  if(nrow(pixel) > window){

    # Loop through moving windows
    for (i in (w_floor + 1):(nrow(pixel) - w_floor)) {
      # Determine observation of interest and window before/after
      center <- pixel[i, ] # Center observation of interest
      pre <- pixel[(i - w_floor):(i - 1), ] # Obs before date of interest
      post <- pixel[(i + 1):(i + w_floor), ] # Obs after date of interest

      # Calculate diffs between the median values before/after central obs
      pre_diff <- center$signal - median(pre$signal)
      post_diff <- median(post$signal) - center$signal
      # If differences before/after central obs > threshold deviations
      # between the median values pre- and post-obs of interest
      if ((abs(pre_diff + post_diff) >= (threshold * sd(c(pre$signal,
          post$signal), na.rm = TRUE))) &
        # and the range of dates is within the timeframe threshold
        (max(post$timestamp) - min(pre$timestamp) <= timeframe)) {
        # If difference before AND after the central obs >= spike amplitude
        if (((pre_diff <= - spike_amp) & (post_diff >= spike_amp)) |
        ((pre_diff >= spike_amp) & (post_diff <= - spike_amp))) {
            # add this observation date to list of spike dates
            spike_dates <- c(spike_dates, pixel[i, timestamp])
        }
    }
  }
  }
    # If there's at least one spike date, add a "spike" column where 1 == spike
    if (length(spike_dates >= 1)) {
        pixel[timestamp %in% spike_dates, "spike"] <- 1 # spike found
    }
    # Return the input data with the added spike column, if any spikes detected
    screen <- pixel[which(spike == 1)]$timestamp
    return(screen)
}


#***********************************
# Grid search function for spikes - spike_grid()
#***********************************
#' Perform a grid search of parameters to test spike filter
#' 
#'@description
#' Takes in a dataframe for a signal and sequences of spike filter
#' parameters and performs a grid search
#' 
#' @param cell The pixel number to run the spike filter on (integer)
#' @param signal_df The dataframe of values for the input signal (data.frame)
#' @param window The length of the moving window (integer)
#' @param threshThe # of std away from median to flag spike (numeric)
#' @param amp Spike amplitude (units of value being compared) (numeric)
#' @param save_plots Option to save time series plots with spikes (boolean)
#' @param ylab The y-axis label used if save_plots=T (character)
#' 
#' @return A list of indices flagged as spikes for the grid search tests
#' @export

# Spike grid search function
spike_grid <- function (cell, signal_df, window=window, thresh=thresh, amp=amp, save_plots=FALSE, ylab) {
    # Subset to individual pixel
    pix <- signal_df[cell,]
    # Run spike filter
    spikes <- c(spike_filter(signal=pix, window=window, threshold=thresh, spike_amp=amp))

    # This will save all plots for all runs, be cautious of setting this to TRUE
    if (save_plots == TRUE) {
      png(file = paste0(plot_directory, paste(cell, ylab, ".png", sep="")), width=1200, height=400)
      plot(hls_sorted, pix, main = paste0("Pixel:", cell), xlab = "Date", ylab=ylab, ylim= c(0,1))
      points(hls_sorted[spikes], pix[spikes], col = "red", pch = 4,lwd = 4, cex=1.5)
      dev.off()
    }
  return(spikes)
}

#***********************************
# Pass/Fail function - run_PassFail()
#***********************************
#' Takes in a set of training pixels and loops through grid search results
#' to pass or fail each grid search test
#' 
#'@description
#' Takes in the previously generated grid search results along with 
#' training data for any number of pixels and returns whether or not
#' each set of grid search runs were able to accurately flag all the
#' necessary spike points.  Pass/Fail is recorded as a binary with 
#' 1 = Pass and 0 = Fail.
#' 
#' @param training_data Training pixels with known spike indices (data.table)
#' @param grid_results Results from spike_grid (list)
#' @param pixel_id A vector of pixel IDs for training pixels (vector)
#' 
#' @return A list of 1s and 0s denoting whether each grid search run passed
#' or failed based on the supplied set of training pixels and a list of number
#' of false positives spikes identified in each test.
#' @export

run_PassFail <- function(training_data, grid_results, pixel_id) {
  # Creates empty data frame to store 1s and 0s
  results <- matrix(nrow=length(grid_results[[1]]), ncol=dim(training_data)[1])
  false_positives <- matrix(nrow=length(grid_results[[1]]), ncol=dim(training_data)[1])

  # Loop through all training pixels
  for (pixel in 1:length(grid_results)){
    # Loop through all grid search tests and store temporary vectors
    for (test in 1:length(grid_results[[1]])) {
      tmp_testing <- c(grid_results[[pixel]][[test]])
      tmp_training <- c(training_data[pixel, 2:dim(training_data)[2]])

      # If training data has no spikes and no spikes were identified - pass
      if ((sum(is.na(tmp_training)) == length(tmp_training)) & (length(tmp_testing) == 0)) {
        results[test, pixel] <- 1 # 1 for pass
        difference <- length(tmp_testing) - length(tmp_training[tmp_training != "NA"])
        false_positives[test, pixel] <- difference

      # If number of indices that intersect equal the number of training - pass
      } else if (length(intersect(tmp_testing, tmp_training)) == length(tmp_training[tmp_training != "NA"]) & 
      length(tmp_training[tmp_training != "NA"]) != 0) {
        results[test, pixel] <- 1 # 1 for pass
        difference <- length(tmp_testing) - length(tmp_training[tmp_training != "NA"])
        false_positives[test, pixel] <- difference
      
      # All other cases - fail
      } else {
        results[test, pixel] <- 0 # 0 for fail
        false_positives[test, pixel] <- NA
      }
    }
  }
  # Format final matrices
  colnames(results) <- pixel_id
  rownames(results) <- 1:length(grid_results[[1]])
  colnames(false_positives) <- pixel_id
  rownames(false_positives) <- 1:length(grid_results[[1]])
  return(list("Pass_Fail"=results, "False_Positives"=false_positives))
}

#***********************************
# Run grid search and return results - runFullGrid()
#***********************************
#' Takes in a signal and returns what parameter sets pass and fail
#' along with the number of false positives
#' 
#'@description
#' Combines the spike_grid() and runPassFail() functions to easily and more
#' efficiently pass in a signal data frame and return a list of parameter tests
#' that passed and ones that failed along with the cumulative number of false
#' positive spikes identified for each pixel and parameter test
#' 
#' #' @param signal_df The data.frame of values for the signal, this is used
#' to calculate the avg number of observations for calculating fp_limit (data.frame)
#' @param training_path Path to .csv file of training pixels (path)
#' @param pixel_id A vector of pixel IDs for training pixels (vector)
#' @param window_params A sequence of paramters to test for window (vector)
#' @param thresh_params A sequence of parameters to test for threshold (vector)
#' @param amp_params A sequence of parameters to test for amplitude (vector)
#' 
#' @return A list with pass/fail results, false positive totals, 
#' and the grid of parameters that were tested
#' @export

runFullGrid <- function(signal_df, training_path, pixel_id, ylab,
                        window_params, thresh_params, amp_params) {
  signal_testing <- list()
  for (i in pixel_id) {
    signal_test <- grid_search(spike_grid, cell=i, signal_df=signal_df, 
                   params=list(window=window_params, thresh=thresh_params, 
                   amp=amp_params), ylab = ylab)
    results <- signal_test$results
    signal_testing <- append(signal_testing, list(results))
  }

  # Save parameter grid
  signal_grid <- signal_test$tests
  signal_grid$test_id <- 1:dim(signal_grid)[1]

  # Load in values for training pixels
  signal_train <- read.csv(training_path)
  # Format into data.table
  signal_training <- as.data.table(signal_train)

  # Run Pass/Fail test
  signal_pass_fail <- run_PassFail(signal_training, signal_testing, pixel_id)
  return(list(signal_pass_fail, signal_grid))
}


#***********************************
# Find best runs function - find_BestRuns()
#***********************************
#' Finds the top optimal paramter combinations
#' 
#'@description
#' Takes in the results generated from run_PassFail() and finds
#' the top parameter runs based on set threshold for success and 
#' false positive rates
#' 
#' @param grid_results Results object generated from run_FullGrid() (list)
#' @param pixel_id A vector of pixel IDs for training pixels (vector)
#' @param success_percent A percent (decimal) denoting what percent of training
#' pixels the run should have successfully captured 
#' i.e. 0.5 means give me all runs that successfully found all outliers in at
#' least 50% of the training pixel test sets (numeric)
#' @param fp_limit A percent (decimal) denoting the percent threshold of 
#' acceptable false positives found within a time series.
#' i.e. 0.25 means give me all runs where the number of false positives is less
#' than 25% of all observations within the time series (numeric)
#' @param signal_df The data.frame of values for the signal, this is used
#' to calculate the avg number of observations for calculating fp_limit
#' 
#' @return A data.table of the best runs based on input thresholds
#' @export

find_BestRuns <- function(grid_results, pixel_id, success_percent, fp_limit, signal_df) {
  # Format results from run_FullGrid() into data.table
  pass_fail_dt <- as.data.table(grid_results[[1]][1])
  # Add parameter test numbers to align with grid
  pass_fail_dt$Test <- 1:dim(pass_fail_dt)[1]
  # Sum the number of successes for each parameter combination
  successes <- rowSums(pass_fail_dt[,1:length(pixel_id)])

  # Format false positives from run_PassFail into data.table
  false_positives <- as.data.table(grid_results[[1]][2])
  # Add parameter test numbers to align with grid
  false_positives$Test <- 1:dim(pass_fail_dt)[1]
  # Set NA values to zero for summing purposes
  false_positives[is.na(false_positives)] <- 0
  # Sum the total number of false positives for each parameter combination
  false_positives <- rowSums(false_positives[,1:length(pixel_id)])

  # Combine all values into one data.frame
  all_values <- cbind(successes, false_positives, 1:dim(pass_fail_dt)[1])
  colnames(all_values) <- c("successes", "false_positives", "test_id")
  # Convert to data.table
  all_values <- as.data.table(all_values)

  # Subest all results based on input thresholds
  # Success threshold
  test_no <- length(pixel_id)
  success_threshold <- test_no * success_percent

  # False positive threshold
  pix_obs <- NULL
  # Calculates average number of observations for each test pixel and averages
  for (pix in 1:length(pixel_id)) {
    pix_sig <- signal_df[pix,]
    no_obs <- length(pix_sig) - sum(is.na(pix_sig))
    pix_obs <- rbind(pix_obs, no_obs)
  }
  avg_obs <- mean(pix_obs)

  # Calculate fp_treshold based on input percent of time series to accept as fp
  fp_threshold <- fp_limit * avg_obs

  # Subset data.table with all values and order
  best_runs <- all_values[all_values$successes >= success_threshold & all_values$false_positives <= fp_threshold]
  best_runs <- best_runs[order(-best_runs$successes, best_runs$false_positives),]

  # Add the parameter values for the top grid tests for easy viewing
  return_runs <- merge(best_runs, grid_results[[2]], by = "test_id")
  # Return data.table of optimal runs according to specified thresholds
  return(return_runs)
}

#*******************************************************************************
# Start running grid search
#*******************************************************************************
# Each signal will be testing with it's own individual grid search to 
# maintain computational efficiency and pinpoint the best parameters
# for each signal

#***********************************
# Clean/filter data
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

# List of pixels to test
test_pixels <- c(47651, 49547, 20709, 50775, 35106, 29041, 64319, 50120, 45079, 22085, 
                 44478, 52469, 54014, 44702, 25157, 31275, 72120, 73761, 75158, 56521)

#*******************************************************************************
# NDVI
#*******************************************************************************
# Run grid search
ndvi_results <- runFullGrid(ndvi_df_sorted, ndvi_train, test_pixels, 
                            ylab="NDVI", window_params=seq(5,13,2), 
                            thresh_params=seq(0.01, 0.1, 0.01), 
                            amp_params=seq(0.01, 0.1, 0.01))

#***********************************
# Find runs with the best combo of params
#***********************************
# Specify success and fp threshold
ndvi_best <- find_BestRuns(ndvi_results, test_pixels, 0.5, 0.7, ndvi_df_sorted)

#***********************************
# Get Plots
#***********************************
# Run spike_grid() for selected parameter set to save the plots
for (i in test_pixels) {
  spike_grid(i, ndvi_df_sorted, window=7, thresh=0.2, amp=0.1, save_plots=T, ylab="NDVI")
}


#*******************************************************************************
# SWIR
#*******************************************************************************
# Run grid search
swir_results <- runFullGrid(swir2_df_sorted, swir_train, test_pixels, 
                            ylab="SWIR", window_params=seq(5,13,2), 
                            thresh_params=seq(0.01, 0.1, 0.01), 
                            amp_params=seq(0.01, 0.1, 0.01))

#***********************************
# Find runs with the best combo of params
#***********************************
# Specify success and fp threshold
swir_best <- find_BestRuns(swir_results, test_pixels, 0.5, .1, swir2_df_sorted)

#***********************************
# Get Plots
#***********************************
# Run spike_grid() for selected parameter set to save the plots
for (i in test_pixels) {
  spike_grid(i, swir2_df_sorted, window=13, thresh=0.08, amp=0.09, save_plots=T, ylab="SWIR")
}


#*******************************************************************************
# High Albedo
#*******************************************************************************
# Run grid search
ha_results <- runFullGrid(high_df_sorted, ha_train, test_pixels, 
                          ylab="High Albedo", window_params=seq(5,13,2), 
                          thresh_params=seq(0.1, 1, 0.1), 
                          amp_params=seq(0.1, 1, 0.1))


#***********************************
# Find runs with the best combo of params
#***********************************
# Specify success and fp threshold
ha_best <- find_BestRuns(ha_results, test_pixels, 0.1, .99, high_df_sorted)

#***********************************
# Get Plots
#***********************************
# Run spike_grid() for selected parameter set to save the plots
for (i in test_pixels) {
  spike_grid(i, high_df_sorted, window=3, thresh=0.00000001, amp=0.0001, save_plots=T, ylab="High Albedo")
}


#*******************************************************************************
# Low Albedo
#*******************************************************************************
# Run grid search
la_results <- runFullGrid(low_df_sorted, la_train, test_pixels, 
                          ylab="Low Albedo", window_params=seq(5,13,2), 
                          thresh_params=seq(0.1, 1, 0.1), 
                          amp_params=seq(0.1, 1, 0.1))

#***********************************
# Find runs with the best combo of params
#***********************************
# Specify success and fp threshold
la_best <- find_BestRuns(la_results, test_pixels, 0.5, .99, low_df_sorted)

#***********************************
# Get Plots
#***********************************
# Run spike_grid() for selected parameter set to save the plots
for (i in test_pixels) {
  spike_grid(i, low_df_sorted, window=7, thresh=0.01, amp=0.1, save_plots=T, ylab="Low Albedo")
}

#*******************************************************************************
# Soil
#*******************************************************************************
# Run grid search
soil_results <- runFullGrid(soil_df_sorted, soil_train, test_pixels, 
                          ylab="Soil", window_params=seq(5,13,2), 
                          thresh_params=seq(0.1, 1, 0.1), 
                          amp_params=seq(0.1, 1, 0.1))

#***********************************
# Find runs with the best combo of params
#***********************************
# Specify success and fp threshold
soil_best <- find_BestRuns(soil_results, test_pixels, 0.5, .99, soil_df_sorted)

#***********************************
# Get Plots
#***********************************
# Run spike_grid() for selected parameter set to save the plots
for (i in test_pixels) {
  spike_grid(i, soil_df_sorted, window=7, thresh=0.01, amp=0.1, save_plots=T, ylab="Soil")
}

#*******************************************************************************
# Vegetation
#*******************************************************************************
# Run grid search
veg_results <- runFullGrid(veg_df_sorted, veg_train, test_pixels, 
                          ylab="Vegetation", window_params=seq(5,13,2), 
                          thresh_params=seq(0.1, 1, 0.1), 
                          amp_params=seq(0.1, 1, 0.1))

#***********************************
# Find runs with the best combo of params
#***********************************
# Specify success and fp threshold
veg_best <- find_BestRuns(veg_results, test_pixels, 0.5, .99, veg_df_sorted)

#***********************************
# Get Plots
#***********************************
# Run spike_grid() for selected parameter set to save the plots
for (i in test_pixels) {
  spike_grid(i, veg_df_sorted, window=7, thresh=0.01, amp=0.1, save_plots=T, ylab="Vegetation")
}



