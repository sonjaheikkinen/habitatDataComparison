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
run_check_parameter_effects <- FALSE
run_check_model_fits <- FALSE

buffer_width <- 100 # how far away from line transect are habitat values read from raster (meters)
number_of_samples <- 500 # number of samples for each chain
thinning_values <- c(1) # with thin value x, only every x:th value from chain is taken as sample
modelfit_folds <- 2 # how many folds to use in cross-validation (1 fold for test, all others are training)
overwrite_modelfits <- FALSE # should existing modelfit values be overwritten

# Old lists for exploratory analysis, TO DO: move outside main?
habitat_data_types <- c("fraction", "landscapemetrics")
species_data_types <- c("fraction", "diversity", "presence")

# DIRECTORY PATHS
current_computer <- "main"
paths_df <- read.csv("filepath_config.txt", sep = ";")
dir_base <- paths_df[paths_df$computer == current_computer & paths_df$path_name == "base", ]$path
dir_data <- file.path(dir_base, "Data")
dir_models <- file.path(dir_base, "Models")
dir_results <- file.path(dir_base, "Results")
dir_fitted <- file.path(dir_models, "Fitted")
dir_modelfits <- file.path(dir_results, "Modelfits")

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
library(cluster) # For silhouette scores
library(corrplot) # For reordering correlation matrices
library(ggrepel) # For better label placement
library(plotly) # For 3D plots
library(gridExtra) # For arranging plots in grid
library(viridis) # For colors

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

# CHECK PARAMETER EFFECTS
if (run_check_parameter_effects) {
    source("check_parameter_effects.R")
}

# CHECK MODEL FITS
if (run_check_model_fits) {
    source("check_model_fits.R")
}

# COMPARE MODEL FITS

# EXPLORE REASONS FOR DIFFERENCES



