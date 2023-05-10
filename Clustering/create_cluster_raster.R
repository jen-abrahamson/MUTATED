#*******************************************************************************
# Project: roboBayes - Post Process Clustering/Filter Dev
# Script Purpose: Developing cluster filter
# Date: 2022-02-22
# Author: Jenna Abrahamson
#*******************************************************************************

# Load all necessary libraries
library(raster)

#***********************************
# Set paths to load results.Rdata and ref.tif
#***********************************
load("./KR_Spike_Results.Rdata")

# Load in ref tif
tmp_s_path <- "./interim/ref.tif"
tmp_s <- raster(tmp_s_path)


#***********************************
# Miscelaneous Functions needed
#***********************************

# Function to calculate amplitude
CalcAmp <- function(sin_mod, cos_mod){
    Amp <- sqrt(sin_mod^2 + cos_mod^2)
    return(Amp)
}

#***********************************
# Generate a matrix of inputs for clustering
#***********************************
# This section puts all model parameters of interest into a matrix

# Function to pull from results() and create matrix
getClusterInfo <- function(results){
    models <- matrix(, nrow = length(results), ncol = 38)
    colnames(models) <- c("NDVI_pre", "NDVI_post", "SWIR_pre", "SWIR_post", "High_pre", 
                        "High_post", "Low_pre", "Low_post", "Soil_pre", "Soil_post", "Veg_pre", "Veg_post", 
                        "NDVI_mod_amp", "NDVI_curr_amp", "SWIR_mod_amp", "SWIR_curr_amp", 
                        "High_mod_amp", "High_curr_amp", "Low_mod_amp", "Low_curr_amp",
                        "Soil_mod_amp", "Soil_curr_amp", "Veg_mod_amp", "Veg_curr_amp", 
                        "NDVI_curr_var", "NDVI_mod_var", "SWIR_curr_var", "SWIR_mod_var", 
                        "High_curr_var", "High_mod_var", "Low_curr_var", "Low_mod_var", 
                        "Soil_curr_var", "Soil_mod_var", "Veg_curr_var", "Veg_mod_var",
                        "CP_date", "CP_detected")
    for(i in 1:length(results)){
    # If models exits (i.e. changepoint exists)
        if(length(results[i[[1]]]$mods) == 0){
            models[i,] <- NA
        }
        if(length(results[i][[1]]$mods) > 0){
        # Populate matrix with model values
            # NDVI model values
            models[i, 1] <- results[i][[1]]$mods[[1]]$B[1,1]
            models[i, 2] <- results[i][[1]]$currentModel$B[1,1]
            # SWIR model values
            models[i, 3] <- results[i][[1]]$mods[[1]]$B[1,2]
            models[i, 4] <- results[i][[1]]$currentModel$B[1,2]
            # High albedo model values
            models[i, 5] <- results[i][[1]]$mods[[1]]$B[1,3]
            models[i, 6] <- results[i][[1]]$currentModel$B[1,3]
            # Low albedo model values
            models[i, 7] <- results[i][[1]]$mods[[1]]$B[1,4]
            models[i, 8] <- results[i][[1]]$currentModel$B[1,4]
            # Soil model values
            models[i, 9] <- results[i][[1]]$mods[[1]]$B[1,5]
            models[i, 10] <- results[i][[1]]$currentModel$B[1,5]
            # Vegetation model values
            models[i, 11] <- results[i][[1]]$mods[[1]]$B[1,6]
            models[i, 12] <- results[i][[1]]$currentModel$B[1,6]

            # NDVI amplitude
            ndvi_mod_sin <- results[i][[1]]$mods[[1]]$B[2,1]
            ndvi_currmod_sin <- results[i][[1]]$currentModel$B[2,1]
            ndvi_mod_cos <- results[i][[1]]$mods[[1]]$B[3,1]
            ndvi_currmod_cos <- results[i][[1]]$currentModel$B[3,1]
            ndvi_amp_mods <- CalcAmp(ndvi_mod_sin, ndvi_mod_cos)
            ndvi_amp_currmods <- CalcAmp(ndvi_currmod_sin, ndvi_currmod_cos)
            models[i, 13] <- ndvi_amp_mods
            models[i, 14] <- ndvi_amp_currmods

            # SWIR amplitude
            swir_mod_sin <- results[i][[1]]$mods[[1]]$B[2,2]
            swir_currmod_sin <- results[i][[1]]$currentModel$B[2,2]
            swir_mod_cos <- results[i][[1]]$mods[[1]]$B[3,2]
            swir_currmod_cos <- results[i][[1]]$currentModel$B[3,2]
            swir_amp_mods <- CalcAmp(swir_mod_sin, swir_mod_cos)
            swir_amp_currmods <- CalcAmp(swir_currmod_sin, swir_currmod_cos)
            models[i, 15] <- swir_amp_mods
            models[i, 16] <- swir_amp_currmods

            # High albedo amplitude
            ha_mod_sin <- results[i][[1]]$mods[[1]]$B[2,3]
            ha_currmod_sin <- results[i][[1]]$currentModel$B[2,3]
            ha_mod_cos <- results[i][[1]]$mods[[1]]$B[3,3]
            ha_currmod_cos <- results[i][[1]]$currentModel$B[3,3]
            ha_amp_mods <- CalcAmp(ha_mod_sin, ha_mod_cos)
            ha_amp_currmods <- CalcAmp(ha_currmod_sin, ha_currmod_cos)
            models[i, 17] <- ha_amp_mods
            models[i, 18] <- ha_amp_currmods

            # Low albedo amplitude
            la_mod_sin <- results[i][[1]]$mods[[1]]$B[2,4]
            la_currmod_sin <- results[i][[1]]$currentModel$B[2,4]
            la_mod_cos <- results[i][[1]]$mods[[1]]$B[3,4]
            la_currmod_cos <- results[i][[1]]$currentModel$B[3,4]
            la_amp_mods <- CalcAmp(la_mod_sin, la_mod_cos)
            la_amp_currmods <- CalcAmp(la_currmod_sin, la_currmod_cos)
            models[i, 19] <- la_amp_mods
            models[i, 20] <- la_amp_currmods

            # Soil amplitude
            soil_mod_sin <- results[i][[1]]$mods[[1]]$B[2,5]
            soil_currmod_sin <- results[i][[1]]$currentModel$B[2,5]
            soil_mod_cos <- results[i][[1]]$mods[[1]]$B[3,5]
            soil_currmod_cos <- results[i][[1]]$currentModel$B[3,5]
            soil_amp_mods <- CalcAmp(soil_mod_sin, soil_mod_cos)
            soil_amp_currmods <- CalcAmp(soil_currmod_sin, soil_currmod_cos)
            models[i, 21] <- soil_amp_mods
            models[i, 22] <- soil_amp_currmods

            # Veg amplitude
            veg_mod_sin <- results[i][[1]]$mods[[1]]$B[2,6]
            veg_currmod_sin <- results[i][[1]]$currentModel$B[2,6]
            veg_mod_cos <- results[i][[1]]$mods[[1]]$B[3,6]
            veg_currmod_cos <- results[i][[1]]$currentModel$B[3,6]
            veg_amp_mods <- CalcAmp(veg_mod_sin, veg_mod_cos)
            veg_amp_currmods <- CalcAmp(veg_currmod_sin, veg_currmod_cos)
            models[i, 23] <- veg_amp_mods
            models[i, 24] <- veg_amp_currmods

            # Variance
            # NDVI variance
            ndvi_currmod_var <- results[i][[1]]$currentModel$V[1,1] 
            ndvi_mod_var <- results[i][[1]]$mods[[1]]$V[1,1]
            models[i, 25] <- ndvi_currmod_var
            models[i, 26] <- ndvi_mod_var

            # SWIR variance
            swir_currmod_var <- results[i][[1]]$currentModel$V[2,2] 
            swir_mod_var <- results[i][[1]]$mods[[1]]$V[2,2]
            models[i, 27] <- swir_currmod_var
            models[i, 28] <- swir_mod_var

            # High albdo variance
            ha_currmod_var <- results[i][[1]]$currentModel$V[3,3] 
            ha_mod_var <- results[i][[1]]$mods[[1]]$V[3,3]
            models[i, 29] <- ha_currmod_var
            models[i, 30] <- ha_mod_var

            # Low albedo variance
            la_currmod_var <- results[i][[1]]$currentModel$V[4,4] 
            la_mod_var <- results[i][[1]]$mods[[1]]$V[4,4]
            models[i, 31] <- la_currmod_var
            models[i, 32] <- la_mod_var

            # Soil variance
            soil_currmod_var <- results[i][[1]]$currentModel$V[5,5] 
            soil_mod_var <- results[i][[1]]$mods[[1]]$V[5,5]
            models[i, 33] <- soil_currmod_var
            models[i, 34] <- soil_mod_var

            # Veg variance
            veg_currmod_var <- results[i][[1]]$currentModel$V[6,6] 
            veg_mod_var <- results[i][[1]]$mods[[1]]$V[6,6]
            models[i, 35] <- veg_currmod_var
            models[i, 36] <- veg_mod_var

            # Changepoint & latentcy dates
            first_cp <- results[i][[1]]$cp_info$cp_dates[1]
            cp_detected <- results[i][[1]]$cp_info$cp_detected[1]
            models[i, 37] <- as.Date(first_cp, origin ="1970-01-01")
            models[i, 38] <- as.Date(cp_detected, origin ="1970-01-01")
        }
    }
    return(models)
}


#***********************************
# User Inputs
#***********************************
# Set these user inputs

# Set this to desired number, then run everything below
number_of_clusters <- 10
# Minimum and maximum size thresholds
min_size <- 8
max_size <- 2000

# Set as output raster filename if desired
file_name <- "./KR_R002_Clusters_Norm"

#***********************************
# Run Clustering - run everything until next subsection
#***********************************
# This will generate a quick plot to verify clusters look correct

# Format matrix and remove NA values
results_models <- getClusterInfo(results)
models_to_cluster <- na.omit(results_models)

# Cluster all values
set.seed(123)
clusters <- kmeans(x=models_to_cluster, number_of_clusters)

# Isolate clusters
cp_clusters <- clusters$cluster

# Sort out index values
index_values <- matrix(, nrow = length(results), ncol = 1)
index_values[,1] <- 1:length(results)

all_values <- cbind(results_models, index_values)
colnames(all_values)[39] <- "Index"
no_na_values <- na.omit(all_values)


# Add cluster numbers to matrix
cluster_to_add <- matrix(, nrow=length(cp_clusters), ncol=1)
cluster_to_add[,1] <- cp_clusters

clustered_data <- cbind(no_na_values, cluster_to_add)
colnames(clustered_data)[40] <- "Cluster_ID"

clusters_and_index <- clustered_data[,39:40]

final_clusters <- merge(all_values, clusters_and_index, by = "Index", all=T)


# Get into raster
xy <- coordinates(tmp_s)
xy <- as.matrix(xy)

plot_clusters <- cbind(final_clusters, xy)

plot(plot_clusters$x, plot_clusters$y, col = plot_clusters$Cluster_ID, pch=20)

#***********************************
# Create cluster plot with size filtering
#***********************************
# This will filter based on size and plot a raster

create_raster <- plot_clusters[,40:42]
clu <- plot_clusters[,40]
create_raster <- cbind(create_raster, clu)
create_raster <- create_raster[,2:4]
raster <- rasterFromXYZ(create_raster, crs=crs(tmp_s))

cluster_raster <- raster$clu
clumped_raster <- clump(cluster_raster, directions=4)

f <- freq(clumped_raster)
f <- as.data.frame(f)
# Separate out rows that are only represented by clumps under 8 pixels or over 2000 pixels
str(which(f$count <= min_size | f$count >= max_size))

# Find corresponding values
str(f$value[which(f$count <= min_size | f$count >= max_size)])

# Put values into a vector of clump ID's to be removed
excludeID <- f$value[which(f$count <= min_size | f$count >= max_size)]

# Set unrealistic sized clusters as NA
clumped_sites <- clumped_raster

cluster_raster[clumped_raster %in% excludeID] <- NA
plot(cluster_raster, col = viridis(number_of_clusters))

#***********************************
# Save Raster
#***********************************

# Saves the output raster plotted in previous step
writeRaster(cluster_raster, file_name, overwrite=T)
