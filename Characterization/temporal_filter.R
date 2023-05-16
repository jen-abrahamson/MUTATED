#*******************************************************************************
# Project: roboBayes pipeline
# Script Purpose: filter by temporal window with ST-DBSCAN
# Date: 2022-12-22
# Author: Jenna Abrahamson
#*******************************************************************************
# Load libraries
library(terra)
library(data.table)
source("functions/st-dbscan.R")
source("functions/post_extract_quantities.R")
source("functions/post_filters.R")

# Load results
results_path <- "./ACC.Rdata"
load(results_path)

# Load ref tif
ref_tif_path <- "./tmp.tif"
ref_tif <- rast(ref_tif_path)$tmp_1

# ~~~~~~~~~~~~~~~~~~~~
# Add spatial info to results
# ~~~~~~~~~~~~~~~~~~~~
# get all cps
cps_all <- get_all_cps(results, dates)

# generate cps_first rast
cps_first <- get_first_cps(results, dates)
cps_first_rast <- ref_tif
values(cps_first_rast) <- NA
values(cps_first_rast) <- c(t(matrix(cps_first, nrow(ref_tif),ncol(ref_tif))))

# filter by spatial polygons (size filter)
polys <- to_polygons(cps_first_rast)
polys <- (small_site_filter(polys, size_thresh = 8000))
cps_filtered <- mask(cps_first_rast, polys)
cps_first_filtered <- values(cps_filtered)

# generate index rast
index_vals <- c(1:length(cps_first_filtered))
index_rast <- ref_tif
values(index_rast) <- NA
values(index_rast) <- index_vals


# convert cps first rast to points
cps_first_pts <- as.points(cps_filtered, values=TRUE, na.rm=FALSE)
cps_first_pts$index <- c(1:length(cps_first_filtered))

# generate index data.frame
index_df <- as.data.frame(values(cps_first_pts))
names(index_df) <- c("cp", "index")
index_df$x <- crds(cps_first_pts, df=TRUE)$x
index_df$y <- crds(cps_first_pts, df=TRUE)$y
for (cp in 1:dim(index_df)[1]) {
    if (index_df[cp,"cp"] %in% NA == TRUE) {
        next
    } else {
        date <- index_df[cp,"cp"]
        index_df[cp, "cp"] <- as.Date(as.character(dates[date]))
    }
}

# ~~~~~~~~~~~~~~~~~~~~
# Run spatiotemporal clustering
# ~~~~~~~~~~~~~~~~~~~~
# run ST-DBSCAN
input_df <- index_df
run_clusters <- stdbscan(input_df$x, input_df$y, input_df$cp, 60, 90, 5)

# append cluster numbers to input_dt
input_df$cluster <- run_clusters[1]$cluster
input_dt <- as.data.table(input_df)
# get back into raster to check results
cluster_data <- as.data.frame(input_dt[input_dt$cluster != 0])

# create column to order by
cluster_vals <- merge(index_df, cluster_data, by="index", all.x=TRUE)
cluster_rast <- ref_tif
values(cluster_rast) <- NA
values(cluster_rast) <- cluster_vals$cluster

cluster_polys <- as.polygons(cluster_rast)
cluster_final <- small_site_filter(cluster_polys, size_thresh=8000)

# convert back to polygons
cps_first_final <- mask(cps_first_rast, cluster_rast)
