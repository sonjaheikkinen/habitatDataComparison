# THIS FILE CONTAINS FUNCTIONS THAT ARE EITHER USED IN THE MAIN SCRIPT, 
# OR USED BY MULTIPLE SUBSCRIPTS

# GLOBAL FUNCTIONS USED BY MULTIPLE SUBSCRIPTS

# MAIN SCRIPT FUNCTIONS
read_climate_data <- function(folder, type) {
    orefix <- ''
    if (type == "temperature") {
        prefix <- "tmon"
    }
    if (type == "rainfall") {
        prefix <- "rrmon"
    }
    raster_files <- list.files(folder, pattern = sprintf("^%s", prefix), full.names = TRUE)
    rasters <- c()
    for (filename in raster_files) {
        rasters <- c(rasters, rast(filename))
    }
    return(rasters)
}
