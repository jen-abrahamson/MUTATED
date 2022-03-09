# Spike Filter Grid Search
*Code used for developing pre-processing filters in roboBayes pipeline*

This repository contains functions used for running a grid search for spike filter parameters.  The grid search is 
completed using training data that was generated via visual interpretation of remotely sensed time series.  The 
training data consists of observation indices that correspond to obvious outlier points in 20 different pixels 
for a ~100 km<sup>2</sup> region in Pyongyang, South Korea.  Training data used to run the supplied script is included as .csv files 
in this repository.  Example time series plots generated using the script are also included here.
