# MUTATED
Code developed for MUTATED research assistantship

## Filters
Code developed for classification of "heavy construction" per SMART's goals as determined by IARPA.  This contains parameter filters targeted at classification of change sites.

## Spike_Filter
Code developed for outlier detection in removal in remote sensing time series.  Able to be used in an online monitoring application, this code contains two different options for spike filters:
* 1. Central Moving Window Spike Filter - Flags spikes based on a set threshold (number of standard deviations away from mean/median) compared to values before **and** after the observation as determined by the user-specified window.
* 2. Lagging Window Spike Filter - Flags spikes based on a set threshold (number of standard deviations away from mean/median_ compared to values **before** the observation as determined by the user-specified window.

## Plotting
Code developed for quick plotting of time series from a harmonized data cube containing Landsat-8 and Sentinel-2 Imagery.
