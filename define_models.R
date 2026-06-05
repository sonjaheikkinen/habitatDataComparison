# SCRIPT FOR DEFINING MODELS FOR HMSC

# LOAD DATA
load(file = file.path(dir_data, "abundance.RData"))
load(file = file.path(dir_data, "occurrence.RData"))
load(file = file.path(dir_data, "spatiotemporal_context.RData"))
load(file = file.path(dir_data, "env_data_natura.RData"))
load(file = file.path(dir_data, "env_data_corine.RData"))
load(file = file.path(dir_data, "trait_data.RData"))
load(file = file.path(dir_data, "phylogeny_data.RData"))
load(file = file.path(dir_data, "fractions_natura.RData"))
load(file = file.path(dir_data, "fractions_corine.RData"))
trait_data <- trait_data[,c("Feeding", "Mig"), drop = FALSE]



# Set abundance to NA where 0 for hurdle model approach
#abundance[abundance == 0] <- NA




# FORMULATE STUDY DESIGN
# Rows: samples, columns: hierarchical categories (year, transect) of each sample
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
#randomlevel_spatial$nfMax <- 3


# Create year-level temporal random effect 
# Each sample belongs to a year. Each year has temporal coordinates year
temporal_coordinates <- data.frame(Year = unique(spatiotemporal_context$Year))
rownames(temporal_coordinates) <- temporal_coordinates$Year
randomlevel_temporal <- HmscRandomLevel(sData = temporal_coordinates)
#randomlevel_temporal$nfMax <- 3


# Create non-spatial random effects for transect and year
randomlevel_transect <- HmscRandomLevel(units = levels(study_design$Transect))
randomlevel_year <- HmscRandomLevel(units = levels(study_design$Year))






# WRITE FORMULAS
# For non-categorical variables:
# If you want the occurrence probability to be able to peak at an intermediate value (bell curve)
# You need to use second order polynomial poly(variable, degree = 2, raw = TRUE)
trait_formula <- as.formula("~Feeding")
formula_natura<- as.formula(sprintf("~Effort:(%s+Temperature+Rainfall+PatchDensity)", 
                                          paste(colnames(fractions_natura), 
                                                collapse = "+")))
formula_corine <- as.formula(sprintf("~Effort:(%s+Temperature+Rainfall+PatchDensity)", 
                                     paste(colnames(fractions_corine), 
                                           collapse = "+")))
forest_formula_natura <- ~Effort:(Luonnonmetsät + Tunturikoivikot +  Lehdot + Tulvametsät +
                                  Temperature + Rainfall + PatchDensity)
forest_formula_corine <- ~Effort:(Havumetsät.kivennäismaalla + 
                                  Sekametsät.kivennäismaalla + Sekametsät.turvemaalla + 
                                  Lehtimetsät.kivennäismaalla +
                                  Havumetsät.kalliomaalla + 
                                  Temperature + Rainfall + PatchDensity)
forest_natura_no_effort_formula <- ~(Luonnonmetsät + Tunturikoivikot +  Lehdot + Tulvametsät +
                                      Temperature + Rainfall + PatchDensity)



# DEEFINE SPIKE-AND-SLAB
# Prior probability for covariate to be included
prob_include_variable <- 0.1
# Species and covariate counts (add 1 for intercept)
number_of_natura_covariates <- length(colnames(fractions_natura)) + 1
number_of_corine_covariates <- length(colnames(fractions_corine)) + 1
number_of_species <- length(colnames(occurrence))
# Select variables separately for each species
variable_selection_params_natura <- list()
for (covariate_number in 1:number_of_natura_covariates) {
    covariate_group <- covariate_number
    species_groups <- 1:number_of_species
    variable_inclusion_probabilities <- rep(prob_include_variable, number_of_species)   
    variable_selection_params_natura[[covariate_number]] <- list(covGroup = covariate_group,
                                                                  spGroup = species_groups,
                                                                  q = variable_inclusion_probabilities)
}
variable_selection_params_corine <- list()
for (covariate_number in 1:number_of_corine_covariates) {
    covariate_group <- covariate_number
    species_groups <- 1:number_of_species
    variable_inclusion_probabilities <- rep(prob_include_variable, number_of_species)   
    variable_selection_params_corine[[covariate_number]] <- list(covGroup = covariate_group,
                                                                 spGroup = species_groups,
                                                                 q = variable_inclusion_probabilities)
}


    



# DEFINE MODELS

probit_natura <- Hmsc(Y = occurrence, 
                      XData = env_data_natura,
                      XFormula = formula_natura,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      #phyloTree = phylogeny_data,
                      distr = "probit",
                      studyDesign = study_design,
                      ranLevels = list("Transect" = randomlevel_spatial,
                                       "Year" = randomlevel_temporal))
probit_corine <- Hmsc(Y = occurrence, 
                      XData = env_data_corine,
                      XFormula = formula_corine,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      #phyloTree = phylogeny_data,
                      distr = "probit",
                      studyDesign = study_design,
                      ranLevels = list("Transect" = randomlevel_spatial,
                                       "Year" = randomlevel_temporal))

probit_natura_forest <- Hmsc(Y = occurrence, 
                      XData = env_data_natura,
                      XFormula = forest_formula_natura,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      #phyloTree = phylogeny_data,
                      distr = "probit",
                      studyDesign = study_design,
                      ranLevels = list("Transect" = randomlevel_spatial,
                                       "Year" = randomlevel_temporal))
forest_natura_no_effort <- Hmsc(Y = occurrence, 
                             XData = env_data_natura,
                             XFormula = forest_natura_no_effort_formula,
                             TrData = trait_data,
                             TrFormula = trait_formula,
                             #phyloTree = phylogeny_data,
                             distr = "probit",
                             studyDesign = study_design,
                             ranLevels = list("Transect" = randomlevel_spatial,
                                              "Year" = randomlevel_temporal))


probit_corine_forest <- Hmsc(Y = occurrence, 
                      XData = env_data_corine,
                      XFormula = forest_formula_corine,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      #phyloTree = phylogeny_data,
                      distr = "probit",
                      studyDesign = study_design,
                      ranLevels = list("Transect" = randomlevel_spatial,
                                       "Year" = randomlevel_temporal))
                      
probit_natura_forest_phyl <- Hmsc(Y = occurrence, 
                           XData = env_data_natura,
                           XFormula = forest_formula_natura,
                           TrData = trait_data,
                           TrFormula = trait_formula,
                           phyloTree = phylogeny_data,
                           distr = "probit",
                           studyDesign = study_design,
                           ranLevels = list("Transect" = randomlevel_spatial,
                                            "Year" = randomlevel_temporal))


probit_corine <- Hmsc(Y = occurrence, 
                      XData = env_data_corine,
                      XFormula = formula_corine,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      #phyloTree = phylogeny_data,
                      distr = "probit",
                      studyDesign = study_design,
                      ranLevels = list("Transect" = randomlevel_spatial,
                                       "Year" = randomlevel_temporal))




lognormal_natura <- Hmsc(Y = abundance, 
                         XData = env_data_natura, 
                         XFormula = formula_natura,
                         TrData = trait_data,
                         TrFormula = trait_formula,
                         #phyloTree = phylogeny_data,
                         distr = "lognormal",
                         studyDesign = study_design,
                         ranLevels = list("Transect" = randomlevel_spatial,
                                          "Year" = randomlevel_temporal))
lognormal_corine <- Hmsc(Y = abundance, 
                         XData = env_data_corine,
                         XFormula = formula_corine,
                         TrData = trait_data,
                         TrFormula = trait_formula,
                         #phyloTree = phylogeny_data,
                         distr = "lognormal",
                         studyDesign = study_design,
                         ranLevels = list("Transect" = randomlevel_spatial,
                                          "Year" = randomlevel_temporal))


#natura_spike_and_slab_q10 <- Hmsc(Y = occurrence, 
#                                XData = env_data_natura,
#                                XFormula = test_formula_natura_frac,
#                                XSelect = variable_selection_params_natura,
#                                distr = "probit",
#                                studyDesign = study_design,
#                                ranLevels = list("Transect" = randomlevel_spatial, 
#                                                 "Year" = randomlevel_temporal))
#corine_spike_and_slab_q10 <- Hmsc(Y = occurrence, 
#                                XData = env_data_corine,
#                                XFormula = test_formula_corine_frac,
#                                XSelect = variable_selection_params_corine,
#                                distr = "probit",
#                                studyDesign = study_design,
#                                ranLevels = list("Transect" = randomlevel_spatial, 
#                                                 "Year" = randomlevel_temporal))






#If needed, test that model works properly:
#print(sprintf("Fitting started %s", date()))
fitted <- sampleMcmc(probit_natura_forest, samples = 3, nChains = 4)
#print(sprintf("Fitting ended %s", date()))


# SAVE MODELS
model_list <- list(probit_natura_forest, probit_corine_forest)
names(model_list) <- c("probit_natura_forest", "probit_corine_forest")
save(model_list, file = file.path(dir_models, "models_unfitted.RData"))




# QUICK TESTING





run_quick_test <- function() {
    
    samples <- 50
    thin <- 100
    n <- 4
    folds <- 2


    quick_test_default <- Hmsc(Y = occurrence, 
                               XData = env_data_natura,
                               XFormula = forest_formula_natura,
                               TrData = trait_data,
                               TrFormula = trait_formula,
                               distr = "probit",
                               studyDesign = study_design,
                               ranLevels = list("Transect" = randomlevel_spatial,
                                                "Year" = randomlevel_temporal))
    quick_test_default$ranLevels$Transect$nfMax = Inf
    quick_test_default$ranLevels$Year$nfMax = Inf
    
    quick_test_4 <- quick_test_default
    quick_test_4$ranLevels$Transect$nfMax = 4
    quick_test_4$ranLevels$Year$nfMax = 4
    
    
    quick_test_3 <- quick_test_default
    quick_test_3$ranLevels$Transect$nfMax = 3
    quick_test_3$ranLevels$Year$nfMax = 3
    
    
    
    print(sprintf("Default started %s", date()))
    fitted_model <- sampleMcmc(quick_test_default, 
                               samples = samples, 
                               thin = thin,
                               transient = ceiling(0.5 * samples * thin),
                               nChains = 4,
                               nParallel = 4, 
                               verbose = 10)
    print(sprintf("Default ended %s", date()))
    print(sprintf("Default cross-validaton started %s", date()))
    model <- fitted_model
    partition <- createPartition(model, nfolds = folds, column = "Transect")
    predicted_values_test_set <- computePredictedValues(model, 
                                                        partition = partition,
                                                        nParallel = 4)
    predictive_power <- evaluateModelFit(hM = model, predY = predicted_values_test_set)
    predictive_power_default <- as.data.frame(predictive_power)
    rownames(predictive_power_default) <- colnames(model$Y)
    print(sprintf("Default cross-validaton ended %s", date()))
    waic_default <- computeWAIC(fitted_model)
    
    
    
    
    print(sprintf("4 started %s", date()))
    fitted_model <- sampleMcmc(quick_test_4, 
                               samples = samples, 
                               thin = thin,
                               transient = ceiling(0.5 * samples * thin),
                               nChains = 4,
                               nParallel = 4, 
                               verbose = 10)
    print(sprintf("4 ended %s", date()))
    print(sprintf("4 cross-validaton started %s", date()))
    model <- fitted_model
    #partition <- createPartition(model, nfolds = folds, column = "Transect")
    predicted_values_test_set <- computePredictedValues(model, 
                                                        partition = partition,
                                                        nParallel = 4)
    predictive_power <- evaluateModelFit(hM = model, predY = predicted_values_test_set)
    predictive_power_4 <- as.data.frame(predictive_power)
    rownames(predictive_power_4) <- colnames(model$Y)
    print(sprintf("4 cross-validaton ended %s", date()))
    waic_4 <- computeWAIC(fitted_model)
    
    
    
    
    print(sprintf("3 started %s", date()))
    fitted_model <- sampleMcmc(quick_test_3, 
                               samples = samples, 
                               thin = thin,
                               transient = ceiling(0.5 * samples * thin),
                               nChains = 4,
                               nParallel = 4, 
                               verbose = 10)
    print(sprintf("3 ended %s", date()))
    print(sprintf("3 cross-validaton started %s", date()))
    model <- fitted_model
    #partition <- createPartition(model, nfolds = folds, column = "Transect")
    predicted_values_test_set <- computePredictedValues(model, 
                                                        partition = partition,
                                                        nParallel = 4)
    predictive_power <- evaluateModelFit(hM = model, predY = predicted_values_test_set)
    predictive_power_3 <- as.data.frame(predictive_power)
    rownames(predictive_power_3) <- colnames(model$Y)
    print(sprintf("3 cross-validaton ended %s", date()))
    waic_3 <- computeWAIC(fitted_model)

    
    
    waic_list <- list()
    rmse_list <- list()
    auc_list <- list()
    tjurr2_list <- list()
    result_list <- list(predictive_power_default, predictive_power_4, predictive_power_3)
    waic_results <- c(waic_default, waic_4, waic_3)
    for (result_number in 1:length(result_list)) {
        
        waic_list[[sprintf("%s", result_number)]] <- waic_results[result_number]
        rmse_list[[sprintf("%s", result_number)]] <- mean(result_list[[result_number]]$RMSE)
        auc_list[[sprintf("%s", result_number)]] <- mean(result_list[[result_number]]$AUC)
        tjurr2_list[[sprintf("%s", result_number)]] <- mean(result_list[[result_number]]$TjurR2)
        
    }
    waic_df <- data.frame(Model = names(waic_list), Value = unlist(waic_list), Metric = "WAIC")
    rmse_df <- data.frame(Model = names(rmse_list), Value = unlist(rmse_list), Metric = "RMSE")
    auc_df <- data.frame(Model = names(auc_list), Value = unlist(auc_list), Metric = "AUC")
    tjurr2_df <- data.frame(Model = names(tjurr2_list), Value = unlist(tjurr2_list), Metric = "TjurR2")
    plot_metric(waic_df)
    plot_metric(rmse_df)
    plot_metric(auc_df)
    plot_metric(tjurr2_df)

}

#run_quick_test()

#posterior <- convertToCodaObject(fitted_model)
#samples <- posterior$Beta
#combined_samples <- as.mcmc.list(samples)
#traceplot(combined_samples[,1])








