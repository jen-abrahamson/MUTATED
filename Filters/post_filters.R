#' CalcAmp - calculate amplitude of fitted model
#' 
#'@description
#' Takes in the sin and cosine model coefficients for the
#' fitted model and calculates the amplitude
#' 
#' @param sin_mod The sin coefficient for the fitted model
#' @param cos_mod The cosine coefficient for the fitted model
#' 
#' @return Amplitude value (numeric)
#' 
#' @export

CalcAmp <- function(sin_mod, cos_mod){
    Amp <- sqrt(sin_mod^2 + cos_mod^2)
    return(Amp)
}


#' param_filter - filter out changepoints not characteristic of construction
#' 
#'@description
#' Takes in the NDVI mean, NDVI amplitude, Veg amplitude, 
#' Low Albedo Variance, and limit thresholds for NDVI amp
#' and High Albedo mean
#' 
#' Based on the logic that the current state of a construction
#' area shows the following attributes:
#' Low NDVI mean and amplitude (shouldn't be any vegetation)
#' Low Vegetation amplitude
#' Low Low Albedo variance
#' 
#' For desert sites:
#' low NDVI amplitude before and after a change points along with
#' high High Albedo mean
#' 
#' 
#' @param results The results from Main.R (Rdata)
#' @param cps_first The output from get_first_cps() 
#' @param filt_ndvi_amp NDVI Amplitude threshold (numeric)
#' @param filt_ndvi_mean NDVI Mean threshold (numeric)
#' @param filt_vege_amp Vege Amplitude threshold (numeric)
#' @param filt_la_var Low Albedo Variance threshold (numeric)
#' @param filt_ndvi_curramp_lim NDVI amp threshold for desert filter (numeric)
#' @param filt_ndvi_amp_lim NDVI pre amp threshold for desert filter (numeric)
#' @param filt_ha_lim High Albedo Mean threshold for desert filter (numeric)
#' @return Screened cps_first object
#' 
#' @export
param_filter <- function(results, cps_first, filt_ndvi_amp, 
                        filt_ndvi_mean, filt_vege_amp, filt_la_var, 
                        filt_ndvi_amp_lim, filt_ndvi_curramp_lim, 
                        filt_ha_lim) {

    # Store correct band order for each signal
    ndvi_sig <- match("ndvi", bands)
    vege_sig <- match("vege", bands)
    ha_sig <- match("high", bands)
    la_sig <- match("low", bands)

    # Loop through and pull all necessary model outputs
    for (pixel in 1:length(results)) {
        if (length(results[pixel][[1]]$mods) == 0) {
            next # Skip pixels with no detected change points
        }
        else {
            # Extract previous NDVI amplitude for pixel
            ndvi_mod_sin <- results[pixel][[1]]$mods[[1]]$B[2, ndvi_sig]
            ndvi_mod_cos <- results[pixel][[1]]$mods[[1]]$B[3, ndvi_sig]
            ndvi_amp_mods <- CalcAmp(ndvi_mod_sin, ndvi_mod_cos)

            # Extract current NDVI amplitude for pixel
            ndvi_currmod_sin <- results[pixel][[1]]$currmod$B[2, ndvi_sig]
            ndvi_currmod_cos <- results[pixel][[1]]$currmod$B[3, ndvi_sig]
            ndvi_amp_currmods <- CalcAmp(ndvi_currmod_sin, ndvi_currmod_cos)

            # Extract NDVI mean for pixel
            ndvi_mean_currmod <- results[pixel][[1]]$currmod$B[1, ndvi_sig]

            # Extract veg amplitude for pixel
            vege_currmod_sin <- results[pixel][[1]]$currmod$B[2, vege_sig]
            vege_currmod_cos <- results[pixel][[1]]$currmod$B[3, vege_sig]
            vege_amp_currmods <- CalcAmp(vege_currmod_sin, vege_currmod_cos)

            # Extract low albedo variance for pixel
            la_var_currmod <- results[pixel][[1]]$currmod$V[4,la_sig]

            # Extract high albedo mean for pixel
            ha_mean_currmod <- results[pixel][[1]]$currmod$B[1, ha_sig]

            # Threshold and screen cps_first
            if (ndvi_amp_currmods > filt_ndvi_amp | ndvi_mean_currmod > filt_ndvi_mean |
                vege_amp_currmods > filt_vege_amp | la_var_currmod > filt_la_var | 
                # Desert filter
                (ndvi_amp_mods < filt_ndvi_amp_lim & 
                ndvi_amp_currmods < filt_ndvi_curramp_lim & 
                ha_mean_currmod < filt_ha_lim)) {
                    cps_first[pixel] <- NA # Set to NA to filter out
            }
        }
    }
    # Return filtered cps_first object
    return(cps_first)
}
