# THIS SCRIPT WAS MADE FOR THESIS PROJECT:
# COMPARING EU HABITAT TYPES AND CORINE LAND COVER AS PREDICTORS FOR BIRD DISTRIBUTION IN NORTHERN LAPLAND
# SONJA HEIKKINEN, 2024, LIFE SCIENCE INFORMATICS

# STUDY QUESTIONS:
# 1. How accurately can bird distributions and diversity be predicted using habitat types?
# 2. How does the accuracy compare to other data (such as Corine)? 
# 3. Does the accuracy differ between species? # If the accuracy (or difference in accuracy between habitat/corine) 
# differs between species, what can explain these differences?

# SETTINGS FOR RUNNING THE SCRIPT
run_preprocess <- FALSE

# DIRECTORY PATHS
current_computer <- "main"
paths_df <- read.csv("filepath_config.txt", sep = ";")
dir_base <- paths_df[paths_df$computer == current_computer & paths_df$path_name == "base", ]$path
dir_data <- file.path(dir_base, "Data")
dir_models <- file.path(dir_base, "Models")
dir_results <- file.path(dir_base, "Results")

# ALL NECESSARY LIBRARIES
library(terra) # For rasters and shapefiles
library(readxl) # For excel files
library(ape) # For phylogenetic trees

# GLOBALLY USED FUNCTIONS
source(file = "common_functions.R")

# SCRIPT STARTS
#################################################################################################################

# PREPROCESS DATA
if (run_preprocess) {
    source(file = "preprocess_data.R")
}

# READ IN DATA
natura_raster <- rast(file.path(dir_data, "natura_2393.tif"))
natura_classification <- read_excel(file.path(dir_data, "Ylalappi_luokitus.xls"))
observations <- read.csv(file.path(dir_data, "observations_preprocessed.csv"),
                         sep = ";")
observations$Transect <- as.character(observations$Transect)
transects_shp <- vect(file.path(dir_data, "transects_preprocessed.shp"))
species_traits <- read_excel(file.path(dir_data, "BirdTraits21112018.xlsx"))
species_alternative_names <- read.csv(file.path(dir_data, "speciesAlternativeNames.txt"), sep = ";")
taxonomy <- read.tree(file.path(dir_data, "tree.txt")) #TO DO: CHOOSE RANDOMLY FROM LIST OF TREES?
temperature_data <- read_climate_data(file.path(dir_data, "ilmastoaineisto"), "temperature")
rainfall_data <- read_climate_data(file.path(dir_data, "ilmastoaineisto"), "rainfall")

# FORMAT DATA

# EXPLORATORY ANALYSIS

# SELECT DATA

# DEFINE MODELS

# FIT MODELS

# CHECK MODEL CONVERGENCE

# CHECK MODEL FITS

# COMPARE MODEL FITS

# EXPLORE REASONS FOR DIFFERENCES




