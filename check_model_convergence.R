# SCRIPT FOR CHECKING MODEL CONVERGENCE

fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

for (i in 1:length(fitted_models)) {
    
    load(fitted_models[i])
    model_name <- strsplit(basename(fitted_models[i]), "\\.")[[1]][1]
    print(model_name)
    
}
