# #################################################################################################
# Title: Post-Processing Clustering Exploration
# Date: 1/17/22
# Author: Jenna Abrahamson
# Description: The first portion goes over loading the data.  After that, theres multiple functions
# defined for cluster exploration.  Demos/testing demonstrations are at the bottom.
# #################################################################################################

# Load necessary libraries
library(terra)
library(raster)
library(viridis)
library(geojsonio)

# Load .RData files

# Load in original data
load("./Region_Data.Rdata")
# Load in results data
load("./Results.Rdata")

# Load in reference tif for site

# File path
tmp_s_path <- "./ref.tif"
# Read in as raster
tmp_s <- raster(tmp_s_path)
# Extract x and y coordinates of each pixel
xy <- xyFromCell(tmp_s, 1:ncell(tmp_s))

# Load in geojson site annotations
annotations <- geojson_read("./Positive_Sites.geojson", what="sp")
annotations <- spTransform(annotations, crs(tmp_s))
# ===================================================================================================
# FUNCTIONS
# ===================================================================================================
# Calculate model differences for all signals and store in one list

change_results <- list(1:length(results))
change_filter <- function(results){
  for(k in 1:length(results)){
    # If models exits (i.e. changepoint exists)
    if(length(results[k][[1]]$mods) > 0){
      # Take the difference between mean model of current model and past model
        ndvi_diff <- (results[k][[1]]$currmod[[1]][1,1] - results[k][[1]]$mods[[1]]$B[1,1])
        swir_diff <- (results[k][[1]]$currmod[[1]][1,2] - results[k][[1]]$mods[[1]]$B[1,2])
        highalb_diff <- (results[k][[1]]$currmod[[1]][1,3] - results[k][[1]]$mods[[1]]$B[1,3])
        lowalb_diff <- (results[k][[1]]$currmod[[1]][1,4] - results[k][[1]]$mods[[1]]$B[1,4])
        soil_diff <- (results[k][[1]]$currmod[[1]][1,5] - results[k][[1]]$mods[[1]]$B[1,5])
        veg_diff <- (results[k][[1]]$currmod[[1]][1,6] - results[k][[1]]$mods[[1]]$B[1,6])
        signals <- list(ndvi_diff, swir_diff, highalb_diff, lowalb_diff, soil_diff, veg_diff)
      # Logic
        #change_results <- append(change_results, list(ndvi_diff, swir_diff, highalb_diff, lowalb_diff, soil_diff, veg_diff))
        #change_results <- append(change_results[k], list(signals))
        change_results[k] <- list(signals)
      }
    else {
        change_results[k] <- NA
     }
  }
  return(change_results)
}

# ===================================================================================================
# Create NDVI Change Magnitude Map

# Function that returns list of change magnitudes for NDVI
ndvi_change_mag <- function(results){
    ndvi_change <- list(1:length(results))
    for(j in 1:length(results)){
        if(length(results[j][[1]]$mods) == 1){
        # Take the difference between mean model of current model and past model
            ndvi_diff <- results[j][[1]]$currentModel[[1]][1,1] - results[j][[1]]$mods[[1]]$B[1,1]
            ndvi_change[j] <- ndvi_diff
            }
        if(length(results[j][[1]]$mods) > 1){
            ndvi_diff <- results[j][[1]]$mods[[2]]$B[1,1] - results[j][[1]]$mods[[1]]$B[1,1]
            ndvi_change[j] <- ndvi_diff
            }
        if(length(results[j][[1]]$mods) == 0){
            ndvi_change[j] <- NA
        }
    }
    return(ndvi_change)
}

# Get into raster
getNDVI_Raster <- function(results, change_mag){
    ndvi_changes <- cbind(xy, change_mag)
    ndvi_df <- as.data.frame(ndvi_changes)
    ndvi_rast <- rasterFromXYZ(ndvi_df)
    plot(tmp_s, col = gray.colors(252), alpha=0.75, legend=F)
    plot(ndvi_rast, col=plasma(252), add=T)
    plot(annotations, lwd=2, add=T)
    return(ndvi_rast)
}

# Get input signal dataframe for clustering
getNDVI_inputdf <- function(results, change_mag){
    ndvi_changes <- cbind(xy, change_mag)
    ndvi_df <- as.data.frame(ndvi_changes)
    ndvi_df$sig_results <- ndvi_df$ndvi_change_mag
    return(ndvi_df)
}

# ===================================================================================================
# Create SWIR2 Change Magnitude Map

# Function that returns list of change magnitudes for SWIR
swir_change_mag <- function(results){
    swir_change <- list(1:length(results))
    for(i in 1:length(results)){
        # If models exits (i.e. changepoint exists)
        if(length(results[i][[1]]$mods) > 0){
            # Take the difference between mean model of current model and past model
            swir_diff <- abs(results[i][[1]]$currmod[[1]][1,2] - results[i][[1]]$mods[[1]]$B[1,2])
            if(abs(swir_diff) <= 0.0005){
                swir_change[i] <- NA
            }
            else{
            swir_change[i] <- swir_diff
        }}
        else {
            swir_change[i] <- NA
        }
    }
    return(swir_change)
}

# Get into raster
getSWIR_Raster <- function(results, change_mag){
    swir_changes <- cbind(xy, change_mag)
    swir_df <- as.data.frame(swir_changes)
    swir_rast <- rasterFromXYZ(swir_df)
    plot(tmp_s, col = gray.colors(252), alpha=0.75, legend=F)
    plot(swir_rast, col=plasma(252), add=T)
    plot(annotations, lwd=2, add=T)
    return(swir_rast)
}

# Get input signal dataframe for clustering
getSWIR_inputdf <- function(results, change_mag){
    swir_changes <- cbind(xy, change_mag)
    swir_df <- as.data.frame(swir_changes)
    swir_df$sig_results <- swir_df$swir_change_mag
    return(swir_df)
}

# ===================================================================================================
# Create High Albedo Change Magnitude Map

# Function that returns list of change magnitudes for high albedo
highalb_change_mag <- function(results){
    highalb_change <- list(1:length(results))
        for(x in 1:length(results)){
        # If models exits (i.e. changepoint exists)
        if(length(results[x][[1]]$mods) > 0){
            # Take the difference between mean model of current model and past model
            highalb_diff <- abs(results[x][[1]]$currmod[[1]][1,3] - results[x][[1]]$mods[[1]]$B[1,3])
            if(abs(highalb_diff) <= 0.00005){
                highalb_change[x] <- NA
            }else{
            highalb_change[x] <- highalb_diff
      } 
        }else{
            highalb_change[x] <- NA
        }
    }
    return(highalb_change)
}

# Get into raster
getHALB_Raster <- function(results, change_mag){
    highalb_changes <- cbind(xy, change_mag)
    highalb_df <- as.data.frame(highalb_changes)
    highalb_rast <- rasterFromXYZ(highalb_df)
    plot(tmp_s, col = gray.colors(252), alpha=0.75, legend=F)
    plot(highalb_rast, col=plasma(252), add=T)
    plot(annotations, lwd=2, add=T)
    return(highalb_rast)
}


# Get input signal dataframe for clustering
getHALB_inputdf <- function(results, change_mag){
    highalb_changes <- cbind(xy, change_mag)
    highalb_df <- as.data.frame(highalb_changes)
    highalb_df$sig_results <- highalb_df$highalb_change_mag
    return(highalb_df)
}

# ===================================================================================================
# Create Low Albedo Change Magnitude Map
# Note: Low albedo seems to work well in Bahrain

lowalb_change_mag <- function(results){
    lowalb_change <- list(1:length(results))
        for(y in 1:length(results)){
        # If models exits (i.e. changepoint exists)
        if(length(results[y][[1]]$mods) > 0){
        # Take the difference between mean model of current model and past model
            lowalb_diff <- abs(results[y][[1]]$currmod[[1]][1,4] - results[y][[1]]$mods[[1]]$B[1,4])
            if(abs(lowalb_diff) < 0.0005){
                lowalb_change[y] <- NA
            }
            else{
            lowalb_change[y] <- lowalb_diff
        }
        }
        else {
            lowalb_change[y] <- NA
        }
    }
    return(lowalb_change)
}

# Get into raster
getLALB_Raster <- function(results, change_mag){
    lowalb_changes <- cbind(xy, change_mag)
    lowalb_df <- as.data.frame(lowalb_changes)
    lowalb_rast <- rasterFromXYZ(lowalb_df)
    plot(tmp_s, col = gray.colors(252), alpha=0.75, legend=F)
    plot(lowalb_rast, col=plasma(252), add=T)
    plot(annotations, lwd=2, add=T)
    return(lowalb_rast)
}

# Get input signal dataframe for clustering
getLALB_inputdf <- function(results, change_mag){
    lowalb_changes <- cbind(xy, change_mag)
    lowalb_df <- as.data.frame(lowalb_changes)
    lowalb_df$sig_results <- lowalb_df$lowalb_change_mag
    return(lowalb_df)
}

# ===================================================================================================
# Create Soil Change Magnitude Map

soil_change_mag <- function(results){
    soil_change <- list(1:length(results))
    for(z in 1:length(results)){
        # If models exits (i.e. changepoint exists)
        if(length(results[z][[1]]$mods) > 0){
        # Take the difference between mean model of current model and past model
            soil_diff <- abs(results[z][[1]]$currmod[[1]][1,5] - results[z][[1]]$mods[[1]]$B[1,5])
            soil_change[z] <- soil_diff
        }
        else {
            soil_change[z] <- NA
        }
    }
    return(soil_change)
}

# Get into raster
soil_results <- soil_change_mag(results)
soil_changes <- cbind(xy, soil_results)
soil_df <- as.data.frame(soil_changes)
soil_rast <- rasterFromXYZ(soil_df)

# Get into raster
getSOIL_Raster <- function(results, change_mag){
    soil_changes <- cbind(xy, change_mag)
    soil_df <- as.data.frame(soil_changes)
    soil_rast <- rasterFromXYZ(soil_df)
    plot(tmp_s, col = gray.colors(252), alpha=0.75, legend=F)
    plot(soil_rast, col=plasma(252), add=T)
    plot(annotations, lwd=2, add=T)
    return(soil_rast)
}

# Get input signal dataframe for clustering
getSOIL_inputdf <- function(results, change_mag){
    soil_changes <- cbind(xy, change_mag)
    soil_df <- as.data.frame(soil_changes)
    soil_df$sig_results <- soil_df$soil_change_mag
    return(soil_df)
}

# ===================================================================================================
# Create Vegetation Change Magnitude Map

veg_change_mag <- function(results){
    veg_change <- list(1:length(results))
    for(v in 1:length(results)){
        # If models exits (i.e. changepoint exists)
        if(length(results[v][[1]]$mods) > 0){
            # Take the difference between mean model of current model and past model
            veg_diff <- abs(results[v][[1]]$currmod[[1]][1,6] - results[v][[1]]$mods[[1]]$B[1,6])
            veg_change[v] <- veg_diff
        }
        else {
            veg_change[v] <- NA
        }
    }
    return(veg_change)
}

# Get into raster
getVEG_Raster <- function(results, change_mag){
    veg_changes <- cbind(xy, change_mag)
    veg_df <- as.data.frame(veg_changes)
    veg_rast <- rasterFromXYZ(veg_df)
    plot(tmp_s, col = gray.colors(252), alpha=0.75, legend=F)
    plot(veg_rast, col=plasma(252), add=T)
    plot(annotations, lwd=2, add=T)
    return(veg_rast)
}

# Get input signal dataframe for clustering
getVEG_inputdf <- function(results, change_mag){
    veg_changes <- cbind(xy, change_mag)
    veg_df <- as.data.frame(veg_changes)
    veg_df$sig_results <- veg_df$veg_change_mag
    return(veg_df)
}

# ===================================================================================================
# Create 6-Signal Change Magnitude Map


signal_change_mag <- function(results){
    signal_change <- list(1:length(results))
    for(s in 1:length(results)){
        # If models exits (i.e. changepoint exists)
        if(length(results[s][[1]]$mods) > 0){
            # Take the difference between mean model of current model and past model
            ndvi_diff <- abs(results[s][[1]]$currmod[[1]][1,1] - results[s][[1]]$mods[[1]]$B[1,1])
            swir_diff <- abs(results[s][[1]]$currmod[[1]][1,2] - results[s][[1]]$mods[[1]]$B[1,2])
            highalb_diff <- abs(results[s][[1]]$currmod[[1]][1,3] - results[s][[1]]$mods[[1]]$B[1,3])
            lowalb_diff <- abs(results[s][[1]]$currmod[[1]][1,4] - results[s][[1]]$mods[[1]]$B[1,4])
            soil_diff <- abs(results[s][[1]]$currmod[[1]][1,5] - results[s][[1]]$mods[[1]]$B[1,5])
            veg_diff <- abs(results[s][[1]]$currmod[[1]][1,6] - results[s][[1]]$mods[[1]]$B[1,6])
            signal_diff <- (ndvi_diff + swir_diff + highalb_diff + lowalb_diff + soil_diff + veg_diff)
            signal_change[s] <- signal_diff
        }
        else {
            signal_change[s] <- NA
        }
    }
    return(signal_change)
}


# Get into raster
getSIG_Raster <- function(results, change_mag){
    sig_changes <- cbind(xy, change_mag)
    sig_df <- as.data.frame(sig_changes)
    sig_rast <- rasterFromXYZ(sig_df)
    plot(tmp_s, col = gray.colors(252), alpha=0.75, legend=F)
    plot(sig_rast, col=plasma(252), add=T)
    plot(annotations, lwd=2, add=T)
    return(sig_rast)
}

# Get input signal dataframe for clustering
getSIG_inputdf <- function(results, change_mag){
    sig_changes <- cbind(xy, change_mag)
    sig_df <- as.data.frame(sig_changes)
    sig_df$sig_results <- sig_df$sig_change_mag
    return(sig_df)
}


# ===================================================================================================
# Average change with NDVI and High Albedo

dual_change_mag <- function(results){
    dual_change <- list(1:length(results))
        for(d in 1:length(results)){
        # If models exits (i.e. changepoint exists)
        if(length(results[d][[1]]$mods) > 0){
        # Take the difference between mean model of current model and past model
            ndvi_diff <- abs(results[d][[1]]$currmod[[1]][1,1] - results[d][[1]]$mods[[1]]$B[1,1])
            highalb_diff <- abs(results[d][[1]]$currmod[[1]][1,3] - results[d][[1]]$mods[[1]]$B[1,3])
            dual_diff <- (ndvi_diff + highalb_diff) / 2
            dual_change[d] <- dual_diff
        }
        else {
            dual_change[d] <- NA
        }
    }
    return(dual_change)
}

# Get into raster
getDUAL_Raster <- function(results, change_mag){
    dual_changes <- cbind(xy, change_mag)
    dual_df <- as.data.frame(dual_changes)
    dual_rast <- rasterFromXYZ(dual_df)
    plot(tmp_s, col = gray.colors(252), alpha=0.75, legend=F)
    plot(dual_rast, col=plasma(252), add=T)
    plot(annotations, lwd=2, add=T)
    return(dual_rast)
}

# Get input signal dataframe for clustering
getDUAL_inputdf <- function(results, change_mag){
    dual_changes <- cbind(xy, change_mag)
    dual_df <- as.data.frame(dual_changes)
    dual_df$sig_results <- dual_df$dual_change_mag
    return(dual_df)
}

# ===================================================================================================

# Load MDST-DBSCAN Function:

########################################################################
# INPUTS :                                                             #
# x     = x-axis coordinate of point data                              #                                  
# y     = y-axis coordinate of point data                              #
# time  = time coordinate of point data                                #
# value = additional attributes of point data for clustering           #
# eps = distance maximum for longitude and latitude                    #
# eps2 =  distance maximum for date                                    #
# eps3 = distance maximum for value                                    #
# minpts = number of points to consider a cluster                      #
########################################################################

mdstdbscan <- function (x, 
                        y, 
                        time,
                        value, 
                        eps, 
                        eps2, 
                        eps3, 
                        minpts) { 

  
  distdata <- cbind.data.frame(x, y)
  time <- time
  value <- value
  
  n <- nrow(distdata)
  
  classn <- cv <- integer(n)
  isseed <- logical(n)
  cn <- integer(1)
  
  for (i in 1:n) {
    unclass <- (1:n)[cv < 1]
    
    ##making distance
    a <- data.frame(x = distdata[i, 1], y = distdata[i, 2])
    fordist <- cbind.data.frame(a, distdata)
    idist <- abs(sqrt((fordist[,1] - fordist[,3])^2 + (fordist[, 2] - 
                                                         fordist[, 4])^2))
    forvaluedist <- cbind.data.frame(value[i], value)
    ivaluedist <- abs(forvaluedist[, 1] - forvaluedist[, 2])
    fortime <- cbind.data.frame(time[i], time)
    itimedist <- abs(fortime[, 1] - fortime[, 2])
    
    if (cv[i] == 0) {
      
      reachables <- intersect(unclass[idist[unclass] <= eps],  
                              unclass[itimedist[unclass] <= eps2])
      reachables <- intersect(reachables, unclass[ivaluedist[unclass] <= eps3])
      if (length(reachables) + classn[i] < minpts)
        cv[i] <- (-1)                    
      else {
        cn <- cn + 1                   
        cv[i] <- cn
        isseed[i] <- TRUE
        reachables <- setdiff(reachables, i)
        unclass <- setdiff(unclass, i)       
        classn[reachables] <- classn[reachables] + 1
        while (length(reachables)) {
          cv[reachables] <- cn           
          ap <- reachables                           
          reachables <- integer()
          
          for (i2 in seq(along = ap)) {
            j <- ap[i2]
            
            ##make distance again when cluster is expanding
            b <- data.frame(x = distdata[j, 1], y = distdata[j, 2])
            jfordist <- cbind.data.frame(b, distdata)
            jdist <- sqrt((jfordist[,1] - jfordist[,3])^2 + 
                            (jfordist[, 2] - jfordist[, 4])^2)
            jforvaluedist <- cbind.data.frame(value[j], value)
            jvaluedist <- abs(jforvaluedist[, 1] - jforvaluedist[, 2])
            jfortime <- cbind.data.frame(time[j], time)
            jtimedist <- abs(jfortime[, 1] - jfortime[, 2])
            
            jreachables <- intersect(unclass[jdist[unclass] <= eps],  
                                     unclass[jtimedist[unclass] <= eps2])
            jreachables <- intersect(jreachables, unclass[jvaluedist[unclass]
                                                          <= eps3])
            
            if (length(jreachables) + classn[j] >= minpts) {
              isseed[j] <- TRUE
              cv[jreachables[cv[jreachables] < 0]] <- cn
              reachables <- union(reachables, jreachables[cv[jreachables] == 0])
            }
            classn[jreachables] <- classn[jreachables] + 1
            unclass <- setdiff(unclass, j)
          }
        }
      }
    }
    if (!length(unclass))
      break
  }
  
  if (any(cv == (-1))) {
    cv[cv == (-1)] <- 0
  }
  result <- list(cluster = cv, eps = eps, 
              eps2 = eps2, eps3 = eps3,
              minpts = minpts, density = classn)
  rm(classn)

  class(result) <- "mdst-dbscan"
  return(result)
}


cluster_models <- function(dates, results, mods_df, eps=90, eps2=12, eps3=0.2, minpts=8){
    hls_sorted <- sort(dates)

    # Return time step index of date that a changepoint was first detected at each pixel
    cps_first <- (sapply(results,function(x){
    if(any(names(x[1]$cp_info)=="cp_dates")){
        df <- cbind(hls_sorted,0)
        if(length(x[1]$cp_info$cp_dates)>=1){
        cp_ind <- which(df[,1]==as.numeric(min(x[1]$cp_info$cp_dates)))[1]
        return(cp_ind)
        }
        else{
        return(NA)
        }
    }
    else{
        return(NA)
    }
    }))

    # Combine signal differences with change date time and include indexes
    signal_cps <- cbind(mods_df, cps_first)
    signal_cps$index <- 1:length(results)

    # Run MDST-DBSCAN Clustering Method
    signal_cps_omit <- na.omit(signal_cps) # Take out NAs - runs very slow if there are NAs

    tmp_x <- unlist(signal_cps_omit$x)
    tmp_y <- unlist(signal_cps_omit$y)
    tmp_time <- unlist(signal_cps_omit$cps_first)
    tmp_sig <- unlist(signal_cps_omit$change_mag)
    clusters <- mdstdbscan(x=tmp_x, y=tmp_y, time=tmp_time, value=tmp_sig, eps=eps, eps2=eps2, eps3=eps3, minpts=minpts)

    cluster_results <- cbind(signal_cps_omit, "cluster" = clusters$cluster)
    raster <- rasterFromXYZ(cluster_results, crs=crs(tmp_s))
    cluster_raster <- raster$cluster
    cluster_raster <- clump(cluster_raster, directions=4)
    f <- freq(cluster_raster)
    f <- as.data.frame(f)

    # Separate out rows that are only represented by clumps under 8 pixels or over 2000 pixels
    str(which(f$count <= 8 | f$count >=2000))

    # Find corresponding values
    str(f$value[which(f$count <= 8 | f$count >= 2000)])

    # Put values into a vector of clump ID's to be removed
    excludeID <- f$value[which(f$count <= 8 | f$count >= 2000)]

    # Set unrealistic sized clusters as NA
    clumped_sites <- cluster_raster
    clumped_sites[cluster_raster %in% excludeID] <- NA


    # Plot new raster
    annotations <- spTransform(annotations, crs(tmp_s))

    plot(tmp_s, col = gray.colors(252), alpha=0.75, legend=F)
    plot(clumped_sites, add=T, legend=F, col=plasma(252))
    plot(annotations, lwd=2, add=T)
    return(clumped_sites)
}



# ===================================================================================================
# Test All Functions

# Find changes in model differences
model_results <- change_filter(results)

# Example for NDVI
ndvi_change_mag <- ndvi_change_mag(results)
ndvi_raster <- getNDVI_Raster(results, ndvi_change_mag)
ndvi_input_df <- getNDVI_inputdf(results, ndvi_change_mag)
ndvi_clusters <- cluster_models(dates, results, ndvi_input_df, eps=90, eps2=10, eps3=0.1, minpts=6)
