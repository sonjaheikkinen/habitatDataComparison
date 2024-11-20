# THIS SCRIPT WAS MADE FOR THESIS PROJECT:
# COMPARING EU HABITAT TYPES AND CORINE LAND COVER AS PREDICTORS FOR BIRD DISTRIBUTION IN NORTHERN LAPLAND
# SONJA HEIKKINEN, 2024, LIFE SCIENCE INFORMATICS

# STUDY QUESTIONS:
# 1. How accurately can bird distributions and diversity be predicted using habitat types?
# 2. How does the accuracy compare to other data (such as Corine)? 
# 3. Does the accuracy differ between species? # If the accuracy (or difference in accuracy between habitat/corine) 
# differs between species, what can explain these differences?

# SETTINGS FOR RUNNING THE SCRIPT

run_preprocess_data <- FALSE
run_format_data <- FALSE
run_exploratory_analysis <- FALSE
run_select_data <- FALSE
run_define_models <- FALSE
run_fit_models <- FALSE
run_check_model_convergence <- FALSE
run_check_model_fits <- FALSE

buffer_width <- 25
habitat_data_types <- c("fraction", "landscapemetrics")
species_data_types <- c("fraction", "diversity", "presence")
number_of_samples <- 250
thinning_values <- c(1, 10, 100)
modelfit_folds <- 2
overwrite_modelfits <- FALSE

# DIRECTORY PATHS
current_computer <- "main"
paths_df <- read.csv("filepath_config.txt", sep = ";")
dir_base <- paths_df[paths_df$computer == current_computer & paths_df$path_name == "base", ]$path
dir_data <- file.path(dir_base, "Data")
dir_models <- file.path(dir_base, "Models")
dir_results <- file.path(dir_base, "Results")
dir_fitted <- file.path(dir_models, "Fitted")
dir_modelfits <- file.path(dir_models, "Modelfits")

# ALL NECESSARY LIBRARIES
library(terra) # For rasters and shapefiles
library(readxl) # For excel files
library(ape) # For phylogenetic trees
library(landscapemetrics) # For calculating landscape metrics 
library(vegan) # For species diversities
library(dendextend) # For clustering
library(pheatmap) # For heatmaps
library(ggplot2) # For plotting
library(Hmsc) # For modeling with Hmsc
library(lme4) # For general linear models
library(vioplot) # For violin plots

# GLOBALLY USED FUNCTIONS
source(file = "common_functions.R")




# SCRIPT STARTS


# Read in rasters
natura_raster <- rast(file.path(dir_data, "natura_2393.tif"))
corine_raster <- rast(file.path(dir_data, "Clc2018_FI20m_2393.tif"))
temperature_data <- read_climate_data(file.path(dir_data, "ilmastoaineisto"), "temperature")
rainfall_data <- read_climate_data(file.path(dir_data, "ilmastoaineisto"), "rainfall")



# PREPROCESS DATA
if (run_preprocess_data) {
    source(file = "preprocess_data.R")
}

# FORMAT DATA
if (run_format_data) {
    source("format_data.R")
}

# EXPLORATORY ANALYSIS
if (run_exploratory_analysis) {
    source("exploratory_analysis.R")
}

# SELECT DATA
if (run_select_data) {
    source("select_data.R")
}

# DEFINE MODELS
if (run_define_models) {
    source("define_models.R")
}

# FIT MODELS
if (run_fit_models) {
    source("fit_models.R")
}

# CHECK MODEL CONVERGENCE
if (run_check_model_convergence) {
    source("check_model_convergence.R")
}

# CHECK MODEL FITS
if (run_check_model_fits) {
    source("check_model_fits.R")
}

# COMPARE MODEL FITS

# EXPLORE REASONS FOR DIFFERENCES



