# MUTATED
<p align="center">
<img width="470" height="200" src="https://github.com/jen-abrahamson/MUTATED/assets/86742376/4f00efb1-50b0-4b2a-a15b-e05f5770b49a">
</p>

Code snippets developed as a part of the *Modeling and Understanding Using Temporal Analysis of Transient Earth Data (MUTATED)* research project.  This research was funded by IARPA as a part of the Space-Based Machine Automated Recognition Technique (SMART) Program.

Many of these code examples build on and work off of the Roboust Online Bayesian Monitoring (roboBayes) algorithm developed by [Wendelberger et al., 2021](https://arxiv.org/abs/2112.12899?context=stat). Some parts of the pipeline code and input data/results are intentionally omitted in this repository due to intellectual property of Accenture Federal Services/IARPA.

## Characterization
Code developed for classification of "heavy construction" per SMART's goals as determined by IARPA.  This contains:
* **1. Heuristic Filters** - Classifies and filters areas of heavy construction based on clustering results and using roboBayes model coefficients.
* **2. Machine Learning Filters** - Classifies and filters heavy construction based on a spectral data cube of inputs. Is developed for use with random forest and extreme gradient boosting models.

## Outliers
Code developed for outlier detection in removal in remote sensing time series.  Able to be used in an online monitoring application, this code contains two different options for spike filters:
* **1. Central Moving Window Spike Filter** - Flags spikes based on a set threshold (number of standard deviations away from mean/median) compared to values before **and** after the observation as determined by the user-specified window. Implementation based on spike filter developed by [Isabella Hinks](https://www.isabellahinks.com/).
* **2. Lagging Window Spike Filter** - Flags spikes based on a set threshold (number of standard deviations away from mean/median) compared to values **before** the observation as determined by the user-specified window.

Also contains an example of a grid search run to determine best spike filter parameters for filtering outliers out of different spectral signals including Linear Spectral Mixing Analysis (LSMA).

## Clustering
Code developed for unsupervised clustering (K-Means) of time series model parameters for detected change sites.  Clustering was used to aid in characterization and classification of construction sites.  Number of inputs dependent on number of spectral signals used for change detection.

## Visualization
Code developed for visualization of MUTATED data inputs, and results. Contains code for: 
* 1. Quick plotting of time series from a harmonized data cube containing Landsat-8 and Sentinel-2 Imagery.
* 2. Creation of animation of detected change sites over time (code for creating .gif)
