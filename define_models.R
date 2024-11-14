# SCRIPT FOR DEFINING MODELS FOR HMSC

# LOAD DATA
load(file = file.path(dir_data, "abundance.RData"))
load(file = file.path(dir_data, "occurrence.RData"))
load(file = file.path(dir_data, "spatiotemporal_context.RData"))
load(file = file.path(dir_data, "env_data_natura.RData"))
load(file = file.path(dir_data, "env_data_corine.RData"))
load(file = file.path(dir_data, "trait_data.RData"))




# FORMULATE STUDY DESIGN
# Rows: samples, columns: hierarchical categories (year, transct) of each sample
# This is pretty much same as spatiotemporal context, but with factors instead
study_design <- as.data.frame(unclass(spatiotemporal_context), stringsAsFactors = TRUE)
study_design$Year <- as.factor(study_design$Year)
study_design$x <- NULL
study_design$y <- NULL





# FORMULATE RANDOM EFFECTS

# Create sample-level random effect (biotic interactions) 
randomlevel_sample <- HmscRandomLevel(units = levels(study_design$Sample))


# Create transect-level spatial random effect
# Each sample belongs to a transect. Each transect has spatial coordinates x and y
# Use sMethod to select an approximation for spatial random level if too slow otherwise
# Use longLat if you don't have metric coordinates
spatial_coordinates <- unique(spatiotemporal_context[, c("Transect", "x", "y")])
rownames(spatial_coordinates) <- spatial_coordinates$Transect
spatial_coordinates$Transect <- NULL
randomlevel_spatial <- HmscRandomLevel(sData = spatial_coordinates)


# Create year-level temporal random effect 
# Each sample belongs to a year. Each year has temporal coordinates year
temporal_coordinates <- data.frame(Year = unique(spatiotemporal_context$Year))
rownames(temporal_coordinates) <- temporal_coordinates$Year
randomlevel_temporal <- HmscRandomLevel(sData = temporal_coordinates)

# Create transect-year-level spatiotemporal random effect
# Each sample belongs to a spatiotemporal class. 
# Each spatiotemporal class has coordinates x, y, year
spatiotemporal_coordinates <- unique(spatiotemporal_context[, c("YearTransect", "x", "y", "Year")])
rownames(spatiotemporal_coordinates) <- spatiotemporal_coordinates$YearTransect
spatiotemporal_coordinates$YearTransect <- NULL
randomlevel_spatiotemporal <- HmscRandomLevel(sData = spatiotemporal_coordinates)






# WRITE FORMULAS
# For non-categorical variables:
# If you want the occurrence probability to be able to peak at an intermediate value (bell curve)
# You need to use second order polynomial poly(variable, degree = 2, raw = TRUE)
env_formula_natura <- as.formula(sprintf("~%s", paste(colnames(env_data_natura), collapse = "+")))
env_formula_corine <- as.formula(sprintf("~%s", paste(colnames(env_data_corine), collapse = "+")))
trait_formula <- as.formula(sprintf("~%s", paste(colnames(trait_data), collapse = "+")))

# DEFINE MODELS
probit_natura <- Hmsc(Y = occurrence, 
                      XData = env_data_natura,
                      XFormula = env_formula_natura,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      phyloTree = taxonomy,
                      distr = "probit",
                      studyDesign = study_design,
                      ranLevels = list("Transect" = randomlevel_spatial, 
                                       "Year" = randomlevel_temporal))

probit_corine <- Hmsc(Y = occurrence, 
                      XData = env_data_corine, 
                      XFormula = env_formula_corine,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      phyloTree = taxonomy,
                      distr = "probit",
                      studyDesign = study_design,
                      ranLevels = list("Transect" = randomlevel_spatial, 
                                       "Year" = randomlevel_temporal))


# If needed, test that model works properly:
# sampleMcmc(probit_natura, samples = 3)


# SAVE MODELS
model_list <- list(probit_natura, probit_corine)
names(model_list) <- c("probit_natura", "probit_corine")
save(model_list, file = file.path(dir_models, "models_unfitted.RData"))


