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
trait_formula <- as.formula("~Feeding")
formula_natura<- as.formula(sprintf("~Effort:(%s+Temperature+Rainfall+PatchDensity)", 
                                          paste(colnames(fractions_natura), 
                                                collapse = "+")))
formula_corine <- as.formula(sprintf("~Effort:(%s+Temperature+Rainfall+PatchDensity)", 
                                     paste(colnames(fractions_corine), 
                                           collapse = "+")))



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

probit_natura_notemporal <- Hmsc(Y = occurrence, 
                      XData = env_data_natura,
                      XFormula = formula_natura,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      #phyloTree = phylogeny_data,
                      distr = "probit",
                      studyDesign = study_design,
                      ranLevels = list("Transect" = randomlevel_spatial))
probit_natura_onlytemporal <- Hmsc(Y = occurrence, 
                      XData = env_data_natura,
                      XFormula = formula_natura,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      #phyloTree = phylogeny_data,
                      distr = "probit",
                      studyDesign = study_design,
                      ranLevels = list("Year" = randomlevel_temporal))
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

probit_natura_no_spatial <- Hmsc(Y = occurrence, 
                      XData = env_data_natura,
                      XFormula = formula_natura,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      #phyloTree = phylogeny_data,
                      distr = "probit",
                      studyDesign = study_design)
probit_corine_no_spatial <- Hmsc(Y = occurrence, 
                      XData = env_data_corine,
                      XFormula = formula_corine,
                      TrData = trait_data,
                      TrFormula = trait_formula,
                      #phyloTree = phylogeny_data,
                      distr = "probit",
                      studyDesign = study_design)


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






# If needed, test that model works properly:
#print(sprintf("Fitting started %s", date()))
#fitted <- sampleMcmc(lognormal_corine, samples = 3, nChains = 4)
#print(sprintf("Fitting ended %s", date()))


# SAVE MODELS
model_list <- list(probit_natura_notemporal, probit_natura_onlytemporal)
names(model_list) <- c("probit_natura_notemporal", "probit_natura_onlytemporal")
save(model_list, file = file.path(dir_models, "models_unfitted.RData"))




# QUICK TESTING





run_quick_test <- function() {
    
    samples <- 50
    thin <- 1
    n <- 4

    env_formula <- ~Effort:(Karut.kirkasvetiset.järvet + Tunturikankaat + Karut.tunturiniityt + 
                                Vaihettumissuot.ja.rantasuot + Aapasuot + Tunturikoivikot + 
                                Palsasuot + Silikaattikalliot + Luonnonmetsät + Lehdot + 
                                Puustoiset.suot + Tulvametsät + Temperature + Rainfall + 
                                PatchDensity)
    env_formula_short <- ~Effort:(Karut.kirkasvetiset.järvet + Tunturikankaat + Karut.tunturiniityt + 
                                      Vaihettumissuot.ja.rantasuot + Aapasuot + Tunturikoivikot + 
                                      Palsasuot + Silikaattikalliot + Luonnonmetsät + Lehdot + 
                                      Puustoiset.suot + Tulvametsät)
    
    quick_test <- Hmsc(Y = occurrence, 
                       XData = env_data_natura,
                       XFormula = env_formula,
                       TrData = trait_data,
                       TrFormula = trait_formula,
                       distr = "probit",
                       studyDesign = study_design,
                       ranLevels = list("Transect" = randomlevel_spatial, 
                                        "Year" = randomlevel_temporal))
    
    quick_test_short <- Hmsc(Y = occurrence, 
                         XData = env_data_natura,
                         XFormula = env_formula_short,
                         TrData = trait_data,
                         TrFormula = trait_formula,
                         distr = "probit",
                         studyDesign = study_design,
                         ranLevels = list("Transect" = randomlevel_spatial, 
                                          "Year" = randomlevel_temporal))
    
    quick_test_notraits <- Hmsc(Y = occurrence, 
                       XData = env_data_natura,
                       XFormula = env_formula,
                       distr = "probit",
                       studyDesign = study_design,
                       ranLevels = list("Transect" = randomlevel_spatial, 
                                        "Year" = randomlevel_temporal))
    
    quick_test_nospatial <- Hmsc(Y = occurrence, 
                       XData = env_data_natura,
                       XFormula = env_formula,
                       TrData = trait_data,
                       TrFormula = trait_formula,
                       distr = "probit",
                       studyDesign = study_design)
    
    quick_test_noyear <- Hmsc(Y = occurrence, 
                       XData = env_data_natura,
                       XFormula = env_formula,
                       TrData = trait_data,
                       TrFormula = trait_formula,
                       distr = "probit",
                       studyDesign = study_design,
                       ranLevels = list("Transect" = randomlevel_spatial))
    
    
    
    print(sprintf("Fitting start %s", date()))
    fitted_model <- sampleMcmc(quick_test, 
                               samples = samples, 
                               thin = thin,
                               transient = 0,
                               nChains = n,
                               nParallel = n)
    print(sprintf("End %s", date()))
    print("")
    
    print(sprintf("Fitting short start %s", date()))
    fitted_model <- sampleMcmc(quick_test_short, 
                                 samples = samples, 
                                 thin = thin,
                                 transient = 0,
                                 nChains = n,
                                 nParallel = n)
    print(sprintf("End %s", date()))
    print("")
    
    print(sprintf("Fitting no traits start %s", date()))
    fitted_model <- sampleMcmc(quick_test_notraits, 
                               samples = samples, 
                               thin = thin,
                               transient = 0,
                               nChains = n,
                               nParallel = n)
    print(sprintf("End %s", date()))
    print("")
    
    print(sprintf("Fitting no year start %s", date()))
    fitted_model <- sampleMcmc(quick_test_noyear, 
                               samples = samples, 
                               thin = thin,
                               transient = 0,
                               nChains = n,
                               nParallel = n)
    print(sprintf("End %s", date()))
    print("")
    
    print(sprintf("Fitting no spatial start %s", date()))
    fitted_model <- sampleMcmc(quick_test_nospatial, 
                               samples = samples, 
                               thin = thin,
                               transient = 0,
                               nChains = n,
                               nParallel = n)
    print(sprintf("End %s", date()))
    print("")
    
    
    
    return(fitted_model)

}

fitted_model <- run_quick_test()

posterior <- convertToCodaObject(fitted_model)
samples <- posterior$Beta
combined_samples <- as.mcmc.list(samples)
traceplot(combined_samples[,1])








