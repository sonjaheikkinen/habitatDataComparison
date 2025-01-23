# SCRIPT FOR DEFINING MODELS FOR HMSC

# LOAD DATA
load(file = file.path(dir_data, "abundance.RData"))
load(file = file.path(dir_data, "occurrence.RData"))
load(file = file.path(dir_data, "spatiotemporal_context.RData"))
load(file = file.path(dir_data, "env_data_natura.RData"))
load(file = file.path(dir_data, "env_data_corine.RData"))
load(file = file.path(dir_data, "trait_data.RData"))
load(file = file.path(dir_data, "fractions_natura.RData"))
load(file = file.path(dir_data, "fractions_corine.RData"))




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

# REMOVED BECAUSE UNSURE IF HMSC SUPPORTS SPATIOTEMPORAL RANDOM EFFECTS
# Create transect-year-level spatiotemporal random effect
# Each sample belongs to a spatiotemporal class. 
# Each spatiotemporal class has coordinates x, y, year
#spatiotemporal_coordinates <- unique(spatiotemporal_context[, c("YearTransect", "x", "y", "Year")])
#rownames(spatiotemporal_coordinates) <- spatiotemporal_coordinates$YearTransect
#spatiotemporal_coordinates$YearTransect <- NULL
#randomlevel_spatiotemporal <- HmscRandomLevel(sData = spatiotemporal_coordinates)






# WRITE FORMULAS
# For non-categorical variables:
# If you want the occurrence probability to be able to peak at an intermediate value (bell curve)
# You need to use second order polynomial poly(variable, degree = 2, raw = TRUE)
env_formula_natura <- as.formula(sprintf("~%s", paste(colnames(env_data_natura), collapse = "+")))
env_formula_corine <- as.formula(sprintf("~%s", paste(colnames(env_data_corine), collapse = "+")))
trait_formula <- as.formula(sprintf("~%s", paste(colnames(trait_data), collapse = "+")))
test_formula_natura_frac <- as.formula(sprintf("~%s", 
                                               paste(colnames(fractions_natura), 
                                                     collapse = "+")))
test_formula_natura_frac_temp <- as.formula(sprintf("~%s", 
                                                    paste(c(colnames(fractions_natura),
                                                            "Temperature"), 
                                                          collapse = "+")))
test_formula_corine_frac <- as.formula(sprintf("~%s", 
                                               paste(colnames(fractions_corine), 
                                                     collapse = "+")))
test_formula_corine_frac_temp <- as.formula(sprintf("~%s", 
                                                    paste(c(colnames(fractions_corine),
                                                            "Temperature"), 
                                                          collapse = "+")))

# DEFINE MODELS

#probit_natura <- Hmsc(Y = occurrence, 
#                      XData = env_data_natura,
#                      XFormula = env_data_natura,
#                      TrData = trait_data,
#                      TrFormula = trait_formula,
#                      phyloTree = taxonomy,
#                      distr = "probit",
#                      studyDesign = study_design,
#                      ranLevels = list("Transect" = randomlevel_spatial,
#                                       "Year" = randomlevel_temporal))

test_probit_natura_frac <- Hmsc(Y = occurrence, 
                                XData = env_data_natura,
                                XFormula = test_formula_natura_frac,
                                distr = "probit",
                                studyDesign = study_design,
                                ranLevels = list("Transect" = randomlevel_spatial, 
                                                 "Year" = randomlevel_temporal))
test_probit_natura_frac_temp <- Hmsc(Y = occurrence, 
                                XData = env_data_natura,
                                XFormula = test_formula_natura_frac_temp,
                                distr = "probit",
                                studyDesign = study_design,
                                ranLevels = list("Transect" = randomlevel_spatial, 
                                                 "Year" = randomlevel_temporal))
test_probit_corine_frac <- Hmsc(Y = occurrence, 
                                XData = env_data_corine,
                                XFormula = test_formula_corine_frac,
                                distr = "probit",
                                studyDesign = study_design,
                                ranLevels = list("Transect" = randomlevel_spatial, 
                                                 "Year" = randomlevel_temporal))
test_probit_corine_frac_temp <- Hmsc(Y = occurrence, 
                                     XData = env_data_corine,
                                     XFormula = test_formula_corine_frac_temp,
                                     distr = "probit",
                                     studyDesign = study_design,
                                     ranLevels = list("Transect" = randomlevel_spatial, 
                                                      "Year" = randomlevel_temporal))
test_probit_natura_clus <- Hmsc(Y = occurrence, 
                                XData = env_data_natura,
                                XFormula = ~Cluster,
                                distr = "probit",
                                studyDesign = study_design,
                                ranLevels = list("Transect" = randomlevel_spatial, 
                                                 "Year" = randomlevel_temporal))



# If needed, test that model works properly:
#sampleMcmc(test_probit_corine_frac, samples = 3)


# SAVE MODELS
model_list <- list(test_probit_natura_frac,
                   test_probit_natura_frac_temp,
                   test_probit_corine_frac,
                   test_probit_corine_frac_temp,)
names(model_list) <- c("test_probit_natura_frac",
                       "test_probit_natura_frac_temp",
                       "test_probit_corine_frac",
                       "test_probit_corine_frac_temp")
save(model_list, file = file.path(dir_models, "models_unfitted_2.RData"))




# QUICK TESTING





quick_test <- Hmsc(Y = occurrence, 
                   XData = env_data_natura,
                   XFormula = ~Cluster,
                   distr = "probit",
                   studyDesign = study_design,
                   ranLevels = list("Transect" = randomlevel_spatial, 
                                    "Year" = randomlevel_temporal))


run_quick_test <- function() {
    
    samples <- 250
    thin <- 1
    n <- 4
    
    print(sprintf("Fitting start %s", date()))
    fitted_model <- sampleMcmc(quick_test, 
                               samples = samples, 
                               thin = thin,
                               transient = ceiling(0.5 * samples * thin),
                               nChains = n,
                               nParallel = n,
                               verbose = 100)
    print(sprintf("End %s", date()))
    print("")
    
    print(sprintf("Cross-validation start %s", date()))
    partition <- createPartition(fitted_model, nfolds = 2) 
    predicted_values <- computePredictedValues(fitted_model, 
                                                   partition = partition, 
                                                   nParallel = n)
    print(sprintf("End %s", date()))
    print("")

    
    
}

#run_quick_test()









