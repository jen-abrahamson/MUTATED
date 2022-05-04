#*******************************************************************************
# Project: roboBayes - Spike Filters
# Script Purpose: Spike Filter functions for both a moving window (central spike
# filter) and a lagging spike filter (recursive spike filter)
# Date: 2022-02-12
# Author: Jenna Abrahamson
#*******************************************************************************
# Load necessary libraries
library(data.table)

#***********************************
# Moving window spike filter (center point) - spike_center()
#***********************************
#' Detect outlier spikes in a time series
#'
#' @description
#' Takes in a single pixel's time series and uses a moving window to flag
#' spikes that are outside of a set standard deviation threshold, or have
#' a difference in spike amplitude that violates a set threshold.
#'
#' @param signal The time series from which to flag spikes (data.frame)
#' @param window The length of the moving window (integer)
#' @param threshold The # of std away from median to flag spike (numeric)
#' @param spike_amp Spike amplitude (units of value being compared) (numeric)
#' @param timeframe Discard if observations span more than x days (integer)
#' @param dates Dates for the time series used to index values (list)
#'
#' @return A vector of pixel indices to be screened from the input time series
#' @export

spike_center <- function(signal, window = 7, threshold = 0.1, spike_amp = 0.2, timeframe = 100, dates = dates) {
  sig_df <- data.table::as.data.table(signal)
  pixel <- cbind(sig_df, "timestamp" = 1:length(dates))
  pixel <- na.omit(pixel)
  w_floor <- floor(window / 2)
  pixel$spike <- 0
  spike_dates <- c()
  if (nrow(pixel) > window) {
    # Loop through moving windows
    print(w_floor)
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
      if ((abs(pre_diff + post_diff) >= (threshold * sd(c(
        pre$signal,
        post$signal
      ), na.rm = TRUE))) &
        # and the range of dates is within the timeframe threshold
        (max(post$timestamp) - min(pre$timestamp) <= timeframe)) {
        # If difference before AND after the central obs >= spike amplitude
        if (((pre_diff <= -spike_amp) & (post_diff >= spike_amp)) |
          ((pre_diff >= spike_amp) & (post_diff <= -spike_amp))) {
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
# Lagging spike filter (compares to previous values) - spike_lag()
#***********************************
#' Detect outlier spikes in a time series
#'
#' @description
#' Takes in a single pixel's time series and uses a moving window to flag
#' spikes that are outside of a set standard deviation threshold, uses a
#' lagging window to compare the next data observation to the window
#' of values that came before that point
#'
#' @param signal The time series from which to flag spikes (data.frame)
#' @param dates Dates for the time series used to index values (list)
#' @param lag The lag of the moving window (# of obs to use) (integer)
#' @param threshold the z-score at which to signal a spike (integer)
#' @param influence the influence of new signals on mean and std (integer)
#'
#' @return A vector of pixel indices to be screened from the input time series
#' @export

# Function for spike detection
spike_lag <- function(signal, lag = 7, threshold = 3, influence = 0.5, dates = dates) {
  na_idx <- which(!is.na(signal))
  y_na <- na.omit(signal)
  signals <- rep(0, length(y_na))
  filteredy_na <- y_na[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y_na[0:lag], na.rm = TRUE)
  stdFilter[lag] <- sd(y_na[0:lag], na.rm = TRUE)
  if (length(y_na) > lag) {
    for (i in (lag + 1):length(y_na)) {
      if (abs(y_na[i] - avgFilter[i - 1]) > threshold * stdFilter[i - 1]) {
        if (y_na[i] > avgFilter[i - 1]) {
          signals[i] <- 1
        } else {
          signals[i] <- -1
        }
        filteredy_na[i] <- influence * y_na[i] + (1 - influence) * filteredy_na[i - 1]
      } else {
        signals[i] <- 0
        filteredy_na[i] <- y_na[i]
      }
      avgFilter[i] <- mean(filteredy_na[(i - lag):i], na.rm = TRUE)
      stdFilter[i] <- sd(filteredy_na[(i - lag):i], na.rm = TRUE)
    }
  }
  signal[na_idx[which(abs(signals) == 1)]] <- NA
  spike_indices <- na_idx[which(abs(signals) == 1)]
  spike_indices <- unname(spike_indices)
  return(spike_indices)
}


#***********************************
# Function to loop through signals and remove spikes
#***********************************
#' Remove identified spikes from time series
#'
#' @description
#' Loops through each pixels time series and sets identified spikes to NA
#'
#' @param signal_df The sorted signal data (data.frame)
#' @param spikes The spike indices (data.table)
#'
#' @return The signal data frame with outlier points set to NA
#' @export


screen_spikes <- function(signal_df, spikes) {
    for (row in 1:nrow(signal_df)) {
        tmp_ind <- c(spikes[[row]])
        tmp_ind <- unlist(tmp_ind)
        signal_df[row, tmp_ind] <- NA
    }
    return(signal_df)
}

