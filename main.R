# THIS SCRIPT WAS MADE FOR THESIS PROJECT:
# COMPARING EU HABITAT TYPES AND CORINE LAND COVER AS PREDICTORS FOR BIRD DISTRIBUTION IN NORTHERN LAPLAND
# SONJA HEIKKINEN, 2024, LIFE SCIENCE INFORMATICS

# STUDY QUESTIONS:
# 1. How accurately can bird distributions and diversity be predicted using habitat types?
# 2. How does the accuracy compare to other data (such as Corine)? 
# 3. Does the accuracy differ between species? # If the accuracy (or difference in accuracy between habitat/corine) 
# differs between species, what can explain these differences?

# SETTINGS FOR RUNNING THE SCRIPT

# DIRECTORY PATHS
current_computer <- "main"
paths_df <- read.csv("filepath_config.txt", sep = ";")
dir_base <- paths_df[paths_df$computer == current_computer & paths_df$path_name == "base", ]$path
dir_data <- file.path(dir_base, "Data")
dir_models <- file.path(dir_base, "Models")
dir_results <- file.path(dir_base, "Results")

# ALL NECESSARY LIBRARIES
library(terra) # For reading rasters

# GLOBALLY USED FUNCTIONS

# SCRIPT STARTS
#################################################################################################################

# PREPROCESS DATA

# READ IN DATA
habitat_types_raster <- rast(file.path(dir_data, "natura_2393.tif"))

# FORMAT DATA

# EXPLORATORY ANALYSIS

# SELECT DATA

# DEFINE MODELS

# FIT MODELS

# CHECK MODEL CONVERGENCE

# CHECK MODEL FITS

# COMPARE MODEL FITS

# EXPLORE REASONS FOR DIFFERENCES




