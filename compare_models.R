# SCRIPT FOR COMPARING MODELFITS

# FUNCTIONS | MODEL COMPARISON

# FUNCTION FOR PLOTTING MODEL PERFORMANCE METRICS
plot_metric <- function(data) {
    metric_name <- data$Metric[1]
    ggplot(data, 
           aes(x = Value, y = reorder(Model, Value))) +
        geom_bar(stat = "identity", fill = "dodgerblue") +
        geom_text(aes(label = round(Value, 3)), hjust = -0.2, size = 4) +
        labs(x = "Value", y = "Model", title = paste("Model Comparison -", metric_name)) +
        theme_minimal()
}

# SCRIPT STARTS

# LOAD IN MODELFITS
modelfit_files <- list.files(dir_modelfits, pattern="*.RData", full.names=TRUE)


# START WITH OVERALL COMPARISONS

# EXTRACT MEAN MODELFIT VALUES FROM EACH FILE
waic_list <- list()
rmse_list <- list()
auc_list <- list()
tjurr2_list <- list()
for (file_number in 1:length(modelfit_files)) {
    
    model_name <- strsplit(basename(modelfit_files[file_number]), "\\.")[[1]][1]
    
    load(modelfit_files[file_number])
    
    waic_list[[model_name]] <- waic
    rmse_list[[model_name]] <- mean(explanatory_power$RMSE)
    auc_list[[model_name]] <- mean(explanatory_power$AUC)
    tjurr2_list[[model_name]] <- mean(explanatory_power$TjurR2)
    
    
}


# CONVERT LISTS TO DATAFRAMES FOR PLOTTING
waic_df <- data.frame(Model = names(waic_list), Value = unlist(waic_list), Metric = "WAIC")
rmse_df <- data.frame(Model = names(rmse_list), Value = unlist(rmse_list), Metric = "RMSE")
auc_df <- data.frame(Model = names(auc_list), Value = unlist(auc_list), Metric = "AUC")
tjurr2_df <- data.frame(Model = names(tjurr2_list), Value = unlist(tjurr2_list), Metric = "TjurR2")


# PLOT METRICS
plot_metric(waic_df)
plot_metric(rmse_df)
plot_metric(auc_df)
plot_metric(tjurr2_df)


rmse_df <- NULL
auc_df<- NULL
tjurr2_df <- NULL

# NEXT COMPARISONS BY SPECIES
for (file_number in 1:length(modelfit_files)) {
    
    model_name <- strsplit(basename(modelfit_files[file_number]), "\\.")[[1]][1]
    
    load(modelfit_files[file_number])
    
    initial_dataframe <- data.frame(species = rownames(explanatory_power))
    
    if (is.null(rmse_df)) {
        rmse_df <- initial_dataframe
    }
    if (is.null(auc_df)) {
        auc_df <- initial_dataframe
    }
    if(is.null(tjurr2_df)) {
        tjurr2_df <- initial_dataframe
    }
    
    rmse_df[,model_name] <- explanatory_power$RMSE
    auc_df[,model_name] <- explanatory_power$AUC
    tjurr2_df[,model_name] <- explanatory_power$TjurR2
    
}

rownames(rmse_df) <- rmse_df$species
rmse_df$species <- NULL
rownames(auc_df) <- auc_df$species
auc_df$species <- NULL
rownames(tjurr2_df) <- tjurr2_df$species
tjurr2_df$species <- NULL

pheatmap(rmse_df)
pheatmap(auc_df)
pheatmap(tjurr2_df)

# FIT LINEAR MODEL TO EXPLAIN FIT USING MODELLING APPROACH
# TO DO: Add species as random effect?
long_rmse <- stack(rmse_df)
long_auc <- stack(auc_df)
long_tjurr2 <- stack(tjurr2_df)

lm_rmse <- lm(values ~ ind, data = long_rmse)
summary(lm_rmse)
lm_auc <- lm(values ~ ind, data = long_auc)
summary(lm_auc)
lm_tjurr2 <- lm(values ~ ind, data = long_tjurr2)
summary(lm_tjurr2)


plot(auc_df$modelfit_probit_corine_thin_300_fitted, 
     auc_df$modelfit_probit_natura_thin_300_fitted)
abline(a = 0, b = 1)


