# Characterization
I was responsbile for developing several different filters and modules related to classifying changes as heavy construction. The roboBayes algorithm functions as a general change detector, therefore to narrow our results and improve F1 scores for heavy construction specifically, I used a variety of heuristics and machine learning models to filter out non-construction related changes.

<p align="center">
  <img width="750" height="400" src="https://github.com/jen-abrahamson/MUTATED/assets/86742376/c4f1e6d4-dc7f-489b-8f38-b52234b0ac01">
  
*The above figure shows the relative change magnitudes for before and after detected changes in key spectral features.*
</p>

## Filtering by Heuristics
### heuristic_filter.R
This script contains functions for performing heuristic filtering of detected change results. This is done as a post-processing step and relies on roboBayes model coefficients and a series of user defined thresholds. Results using this filter are shown below.

<p align="center">
  <img width="800" height="300" src="https://github.com/jen-abrahamson/MUTATED/assets/86742376/27efb250-1207-42dc-8815-8ded6b4474ec">
  
*The above figure shows results before and after heuristic-based filtering.*
</p>

## Filtering via Spatiotemporal Clustering
### temporal_filter.R, st-dbscan.R
These scripts contain code to filter results based on a spatial and temporal moving window. Only results within a set spatial and temporal distance will be kept.

<p align="center">
  <img width="875" height="450" src="https://github.com/jen-abrahamson/MUTATED/assets/86742376/7fccfa52-794d-44be-a09d-e904ca80cb9b">
  
*The above figure shows results using ST-DBSCAN to perform spatiotemporal clustering.*
</p>
