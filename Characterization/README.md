# Characterization
I was responsbile for developing several different filters and modules related to classifying changes as heavy construction. The roboBayes algorithm functions as a general change detector, therefore to narrow our results and improve F1 scores for heavy construction specifically, I used a variety of heuristics and machine learning models to filter out non-construction related changes.

<p align="center">
  <img width="900" height="500" src="https://github.com/jen-abrahamson/MUTATED/assets/86742376/c4f1e6d4-dc7f-489b-8f38-b52234b0ac01">
  
*The above figure shows the relative change magnitudes for before and after detected changes in key spectral features.*
</p>

## heuristic_filter.R
This script contains functions for performing heuristic filtering of detected change results. This is done as a post-processing step and relies on roboBayes model coefficients and a series of user defined thresholds.
