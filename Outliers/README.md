# Outlier Detection
<p align="center">
  <img width="1100" height="350" src="https://github.com/jen-abrahamson/MUTATED/assets/86742376/0d124bb4-64c4-4a07-9298-39dfa8ad3fad)">
  
  *Above is a SWIR2 time series with an outlier removed where the red "x" is displayed.*
</p>

## Moving Window Spike Filters
### spike_filters.R, spike_gridsearch.R
These scripts contain functions to run both a lagging and center moving window spike filter for online change detection algorithms. There is also code showing an example of a grid search conducted for different remote sensing signals in order to determine the optimal parameters for use with the spike filters. The center spike filter is adapted from a spike filter developed by Isabella Hinks.
