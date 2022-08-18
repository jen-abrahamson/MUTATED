# MUTATED
Code snippets developed as a part of the Modeling and Understanding Using Temporal Analysis of Transient Earth Data (MUTATED) research team.  This research was funded by IARPA as a part of the Space-Based Machine Automated Recognition Technique (SMART) Program.

## Filters
Code developed for classification of "heavy construction" per SMART's goals as determined by IARPA.  This contains heuristic parameter filters I developed for targeted classification of construction change sites.

## Spike_Filter
Code developed for outlier detection in removal in remote sensing time series.  Able to be used in an online monitoring application, this code contains two different options for spike filters:
* 1. Central Moving Window Spike Filter - Flags spikes based on a set threshold (number of standard deviations away from mean/median) compared to values before **and** after the observation as determined by the user-specified window.
* 2. Lagging Window Spike Filter - Flags spikes based on a set threshold (number of standard deviations away from mean/median) compared to values **before** the observation as determined by the user-specified window.

## Plotting
Code developed for quick plotting of time series from a harmonized data cube containing Landsat-8 and Sentinel-2 Imagery.
