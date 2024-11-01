# THIS SCRIPTS FORMATS THE DATA IN SUITABLE FORM FOR HMSC ANALYSES

# FUNCTIONS
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

# SCRIPT STARTS
###################################################################################################################

# READ IN DATA
natura_raster <- rast(file.path(dir_data, "natura_2393.tif"))
natura_classification <- read_excel(file.path(dir_data, "Ylalappi_luokitus.xls"))
observations <- read.csv(file.path(dir_data, "observations_preprocessed.csv"), sep = ";")
transects_shp <- vect(file.path(dir_data, "transects_preprocessed.shp"))
species_traits <- read_excel(file.path(dir_data, "BirdTraits21112018.xlsx"))
species_alternative_names <- read.csv(file.path(dir_data, "speciesAlternativeNames.txt"), sep = ";")
taxonomy <- read.tree(file.path(dir_data, "tree.txt")) #TO DO: CHOOSE RANDOMLY FROM LIST OF TREES?
temperature_data <- read_climate_data(file.path(dir_data, "ilmastoaineisto"), "temperature")
rainfall_data <- read_climate_data(file.path(dir_data, "ilmastoaineisto"), "rainfall")


# GENERAL REFORMATTING OF RAW DATA
names_natura <- data.frame(value = natura_classification$Value, name = natura_classification$NaturaTyyppi)
observations$Transect <- as.character(observations$Transect)



# FORMAT ENVIRONMENTAL DATA X
# Rows: samples, columns: env data of each sample

# Get fractions of each habitat type on each transect
fractions_natura <- get_transect_habitat_data(buffer_width, 
                                              natura_raster, 
                                              transects_shp,
                                              names_natura,
                                              "fraction")
fractions_natura <- fractions_natura[,!is.na(colnames(fractions_natura))]
fractions_natura <- fractions_natura[,sapply(fractions_natura, function(x) var(x, na.rm = TRUE) != 0)]

# FORMAT COMMUNITY DATA Y

# FORMAT TRAIT DATA T

# FORMAT PHYLOGENETIC DATA C

# FORMAT SPATIOTEMPORAL CONTEXT S

# SAVE DATA
# Save also raw data as .RData