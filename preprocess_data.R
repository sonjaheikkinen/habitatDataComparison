# THIS SCRIPT PREPROCESSES DATA: 
# The data is limited to only the study area
# The data is cleaned to remove problematic observations


# FUNCTIONS FOR PREPROCESSING
get_transects_that_overlap_raster <- function(transects, raster) {
    raster_values_for_transects <- extract(raster, transects)
    transects_with_no_NA_values <- c()
    for (transect_id in unique(raster_values_for_transects$ID)) {
        values_for_transect <- raster_values_for_transects[raster_values_for_transects$ID == transect_id,]$natura_2393
        if (sum(is.na(values_for_transect)) < 0.1 * length(values_for_transect)) {
            transects_with_no_NA_values <- c(transects_with_no_NA_values, transect_id)
        }
    }
    filtered_transects <- subset(transects,
                                 transects$ID %in% transects_with_no_NA_values)
    
    return(filtered_transects)
}

# SCRIPT STARTS
###################################################################################################################

# READ DATA
transects_shp <- vect(file.path(dir_data, "bird_transects.shp"))
observations <- read.csv(file.path(dir_data, "birdtransects_toRecbase.csv"))

# LIMIT DATA TO STUDY AREA
transects_shp$ID <- 1:nrow(transects_shp) # Adding numerical id, because existing chr identifier caused issues
transects_that_overlap_natura_raster <- get_transects_that_overlap_raster(transects_shp, natura_raster)

# CLEAN DATA

# SAVE PREPROCESSED DATA