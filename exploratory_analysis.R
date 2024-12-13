# HELPER FUNCTIONS

get_clusters <- function(data, cluster_amount) {
    data$Cluster <- NULL
    clustering <- hclust(dist(data), method = "complete")
    cluster_groups <- cutree(tree = as.dendrogram(clustering), k = cluster_amount)
    return(cluster_groups)
}

plot_cluster_elbow <- function(data, data_name) {
    max_clusters <- 15
    within_cluster_sum_of_squares <- sapply(1:max_clusters,
                                            function(cluster_amount) {
                                                clusters <- kmeans(data, 
                                                                   cluster_amount,
                                                                   nstart = 50,
                                                                   iter.max = 15,
                                                                   )
                                                return(clusters$tot.withinss)
                                            })
    plot(1:max_clusters,
         within_cluster_sum_of_squares,
         type = "b",
         pch = 19,
         frame = FALSE,
         xlab = "Number of clusters",
         ylab = "Total within-clusters sum of squares",
         main = sprintf("Elbow curve for %s", data_name))
}


scale_to_range <- function(values, new_min, new_max) {
    old_min <- min(values)
    old_max <- max(values)
    scaled_values <- new_min + (values - old_min) * (new_max - new_min) / (old_max - old_min)
    return(scaled_values)
}






# EXPLORATORY ANALYSIS FUNCTIONS

explore_bird_data <- function(occurrence, abundance, spatiotemporal_context, coordinates, type) {
    
    pdf(file.path(dir_results, sprintf("exploratory_analysis_observations_%s.pdf", type)))
    
    
    # Total abundance in samples
    sample_abundances <- rowSums(abundance, na.rm = TRUE)
    hist(sample_abundances,
         main = "Histogram of total abundance in samples",
         xlab = sprintf("Abundance, min %s, max %s", min(sample_abundances), max(sample_abundances)),
         ylab = sprintf("Frequency, %s total samples", nrow(abundance)))
    
    # Total species richness in samples
    sample_occurrences <- rowSums(occurrence, na.rm = TRUE)
    hist(sample_occurrences,
         main = "Histogram of species richness in samples",
         xlab = sprintf("Occurrence, min %s, max %s", min(sample_occurrences), max(sample_occurrences)),
         ylab = sprintf("Frequency, total samples: %s, total richness. %s", 
                        nrow(occurrence),
                        ncol(occurrence)))
    
    
    # Adjust margin to fit species names to the left
    old_par <- par(no.readonly = TRUE) 
    par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))
    
    # How many samples was each species found in?
    species_sample_occurrences <- sort(colSums(occurrence))
    barplot(species_sample_occurrences, 
            horiz = TRUE, 
            las = 1, 
            main = "How many samples was each species observed in",
            xlab = "Samples",
            cex.names = 0.7)
    
    # Most common species by occurrence
    common_species_occurrences <- tail(species_sample_occurrences, 10)
    barplot(common_species_occurrences, 
            horiz = TRUE, 
            las = 1, 
            main = "How many samples were common species observed in",
            xlab = "Samples",
            cex.names = 0.7)
    
    # Most rare species by occurrence
    rare_species_occurrences <- species_sample_occurrences[species_sample_occurrences < 5]
    barplot(rare_species_occurrences, 
            horiz = TRUE, 
            las = 1, 
            main = "How many samples were rare species observed in",
            xlab = "Samples",
            cex.names = 0.7)
    
    # Use old margins again
    par(old_par)
    
    
    # Spatial and temporal patterns
    
    # Calculate sample level mean abundances across species
    sample_total_abundances <- matrix(0,
                                      nrow = length(unique(spatiotemporal_context$Year)),
                                      ncol = length(unique(spatiotemporal_context$Transect)))
    rownames(sample_total_abundances) <- rev(sort(unique(spatiotemporal_context$Year)))
    colnames(sample_total_abundances) <- unique(spatiotemporal_context$Transect)
    for (row in 1:nrow(spatiotemporal_context)) {
        year <- spatiotemporal_context[row, ]$Year
        transect <- spatiotemporal_context[row, ]$Transect
        sample <- sprintf("%s%s", year, transect)
        year_index <- match(year, rownames(sample_total_abundances))
        transect_index <- match(transect, colnames(sample_total_abundances))
        row_abundances <- abundance[row, ]
        row_abundances[is.na(row_abundances)] <- 0 
        sample_total_abundances[year_index, transect_index] <- sum(row_abundances)
    }
    pheatmap(sample_total_abundances,
             cluster_rows = FALSE,
             main = "Total abundance in each sample")
    
    
    # Spatial variation in transect total abundance 
    # Sample total abundances: sum of species abundances in each sample
    # Transect mean total abundances: mean of sample total abundances from that transect
    transect_mean_total_abundances <- apply(sample_total_abundances, 
                                            2, 
                                            function(col) {
                                                non_zero_values <- col[col != 0]
                                                if (length(non_zero_values) > 0) {
                                                    mean(non_zero_values)
                                                } else {
                                                    NA  # Return NA if there are no non-zero values
                                                }
                                            })
    spatial_plot_data <- data.frame(Transect = names(transect_mean_total_abundances),
                                    Abundance = transect_mean_total_abundances)
    spatial_plot_data$x <- coordinates[spatial_plot_data$Transect, ]$x
    spatial_plot_data$y <- coordinates[spatial_plot_data$Transect, ]$y
    plot(spatial_plot_data$x, spatial_plot_data$y, 
         cex = spatial_plot_data$Abundance / max(spatial_plot_data$Abundance) * 1.5, 
         pch = 21, 
         col = "black", # Dot appearance
         xlab = "X", 
         ylab = "Y", 
         main = "Transect mean total abundances")
    
    
    # Temporal variation in total abundance
    # Sample total abundance: sum of species abundances in each sample
    # Year mean total abundances: mean of sample total abundancies from that year
    year_mean_total_abundances <- apply(sample_total_abundances, 
                                   1, 
                                   function(row) {
                                       non_zero_values <- row[row != 0]
                                       if (length(non_zero_values) > 0) {
                                           mean(non_zero_values)
                                       } else {
                                           NA  # Return NA if there are no non-zero values
                                       }
                                   })
    names(year_mean_total_abundances) <- rownames(sample_total_abundances)
    barplot(year_mean_total_abundances,
            horiz = TRUE,
            las = 1,
            main = "Year mean total abundances")
    
    
    
    # Spatial variation in species composition
    
    # Create spatial occurrence and abundance data
    transect_occurrence_data_list <- list()
    transect_abundance_data_list <- list()
    for (transect in unique(spatiotemporal_context$Transect)) {
        transect_indices <- match(transect, spatiotemporal_context$Transect)
        transect_occurrences <- occurrence[transect_indices, , drop = FALSE]
        transect_abundances <- abundance[transect_indices, , drop = FALSE]
        species_occurrences_in_transect <- apply(transect_occurrences,
                                                 2,
                                                 function(col) {
                                                     any(col == 1) * 1
                                                 })
        species_abundances_in_transect <- apply(transect_abundances,
                                                2,
                                                function(col) {
                                                    sum(col)
                                                })
        transect_occurrence_data_list[transect] <- list(species_occurrences_in_transect)
        transect_abundance_data_list[transect] <- list(species_abundances_in_transect)
    }
    transect_occurrence_data <- as.data.frame(do.call(rbind, transect_occurrence_data_list))
    transect_abundance_data <- as.data.frame(do.call(rbind, transect_abundance_data_list))
    transect_occurrence_clusters <- get_clusters(transect_occurrence_data, 5)
    transect_species_richness <- rowSums(transect_occurrence_data)
    
    
    # Plot species clusters over transects
    plot(coordinates[rownames(transect_occurrence_data), ]$x,
         coordinates[rownames(transect_occurrence_data), ]$y, 
         cex = 1, 
         pch = 21,
         col = transect_occurrence_clusters,
         xlab = "X", 
         ylab = "Y", 
         main = "Transect species cluster")
    
    
    # Plot species richness over transects
    plot(coordinates[rownames(transect_occurrence_data), ]$x,
         coordinates[rownames(transect_occurrence_data), ]$y, 
         cex = transect_species_richness * 0.05, 
         pch = 21,
         col = "Black",
         xlab = "X",
         ylab = "Y", 
         main = "Transect species richness")
    
    
    # Compare mean species diversity over transects
    sample_simpson_diversities <- diversity(abundance, index = "simpson")
    sample_shannon_diversities <- diversity(abundance, index = "shannon")
    transect_mean_simpson_diversities <- c()
    transect_mean_shannon_diversities <- c()
    for (transect in unique(spatiotemporal_context)$Transect) {
        transect_indices <- match(transect, spatiotemporal_context$Transect)
        transect_mean_simpson_diversities[transect] <- mean(sample_simpson_diversities[transect_indices])
        transect_mean_shannon_diversities[transect] <- mean(sample_shannon_diversities[transect_indices])
    }
    
    # Smaller dot means bigger diversity
    scaled_simpson_diversities <- scale(transect_mean_simpson_diversities)
    plot(coordinates[names(transect_mean_simpson_diversities), ]$x,
         coordinates[names(transect_mean_simpson_diversities), ]$y, 
         cex =  (scaled_simpson_diversities + abs(min(scaled_simpson_diversities))) * 0.3,
         pch = 21,
         col = "Black",
         xlab = "X", 
         ylab = "Y", 
         main = "Transect mean Simpson's diversities")
    
    # Bigger dot means bigger diversity
    scaled_shannon_diversities <- scale(transect_mean_shannon_diversities)
    plot(coordinates[names(transect_mean_shannon_diversities), ]$x,
         coordinates[names(transect_mean_shannon_diversities), ]$y, 
         cex =  (scaled_shannon_diversities + abs(min(scaled_shannon_diversities))) * 0.3, 
         pch = 21,
         col = "Black",
         xlab = "X", 
         ylab = "Y", 
         main = "Transect mean Shannon's diversities")
    
    diversity_correlation <- cor(transect_mean_simpson_diversities, 
                                 transect_mean_shannon_diversities)
    plot(coordinates[names(transect_mean_shannon_diversities), ]$x,
         coordinates[names(transect_mean_shannon_diversities), ]$y, 
         cex =  (scaled_shannon_diversities + abs(min(scaled_shannon_diversities))) * 0.3, 
         col = (scaled_simpson_diversities + abs(min(scaled_simpson_diversities))),
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = sprintf("Both diversities, correlation %s", diversity_correlation))
    
    
    richness_diversity_correlation <- cor(transect_species_richness, 
                                          transect_mean_simpson_diversities)
    plot(coordinates[names(transect_mean_shannon_diversities), ]$x,
         coordinates[names(transect_mean_shannon_diversities), ]$y, 
         cex =  transect_species_richness * 0.05, 
         col = (scaled_simpson_diversities + abs(min(scaled_simpson_diversities))),
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = sprintf("Richness (size) and Simpson's (color), correlation %s", 
                        richness_diversity_correlation))
    
    # Are there truly different transects based on species composition?
    plot_cluster_elbow(transect_occurrence_data, "transect species occurrences")
    
    # Are there truly different transects based on transect abundance
    plot_cluster_elbow(transect_abundance_data, "transect species abundances")
    
    
    dev.off()
    
}


explore_habitat_data <- function(natura, 
                                 corine, 
                                 fractions_natura,
                                 fractions_corine,
                                 spatiotemporal_info, 
                                 coordinates, 
                                 type) {
    
    
    # Remove entries before 2006 becaues there is no temperature_data
    selected_years <- spatiotemporal_context$Year >= 2006
    natura <- natura[selected_years,]
    corine <- corine[selected_years,]
    spatiotemporal_context <- spatiotemporal_info[selected_years,]
    
    print(str(natura))
    print(str(spatiotemporal_context))
    
    
    # Create aggregated data where each transect is only once
    data_to_aggregate <- data.frame(
        Transect = spatiotemporal_context$Transect,
        natura
    )
    transect_data_natura <- aggregate(. ~ Transect, data = data_to_aggregate, FUN = mean)
    data_to_aggregate <- data.frame(
        Transect = spatiotemporal_context$Transect,
        corine
    )
    transect_data_corine <- aggregate(. ~ Transect, data = data_to_aggregate, FUN = mean)
    
    
    
    pdf(file.path(dir_results, sprintf("exploratory_analysis_habitats_%s.pdf", type)))

    # Plot histograms of natura and corine types in transects
    # Also plot the histogram of just those values that aren't zeros
    for (column in colnames(fractions_natura)) {
        par(mfrow = c(2, 1))
        column_data <- fractions_natura[, column]
        column_data_without_zeros <- column_data[column_data != 0]
        zeros <- sum(column_data == 0)
        hist(fractions_natura[,column], 
             probability = TRUE, 
             col = "lightgray",
             xlab = sprintf("%s", column),
             main = sprintf("Natura: Histogram of %s in transects, zeros: %s/%s = %s", 
                            column, 
                            zeros, 
                            length(column_data), 
                            round(zeros / length(column_data), 2)))
        if (length(column_data_without_zeros) > 0) {
            hist(column_data_without_zeros, 
                 probability = TRUE, 
                 col = "lightgray",
                 xlab = sprintf("%s", column),
                 main = sprintf("Natura: Histogram of %s without zeros", column))
        }
    }
    for (column in colnames(fractions_corine)) {
        par(mfrow = c(2, 1))
        column_data <- fractions_corine[, column]
        column_data_without_zeros <- column_data[column_data != 0]
        zeros <- sum(column_data == 0)
        hist(fractions_corine[,column], 
             probability = TRUE, 
             col = "lightgray",
             xlab = sprintf("%s", column),
             main = sprintf("Corine: Histogram of %s in transects, zeros: %s/%s = %s", 
                            column, 
                            zeros, 
                            length(column_data), 
                            round(zeros / length(column_data), 2)))
        if (length(column_data_without_zeros) > 0) {
            hist(column_data_without_zeros, 
                 probability = TRUE, 
                 col = "lightgray",
                 xlab = sprintf("%s", column),
                 main = sprintf("Corine: Histogram of %s without zeros", column))
        }
    }
    
    par(mfrow=c(1,1))
    
    hist(natura$Temperature, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Temperature (april + may)",
         main = "Histogram of transect spring temperatures in samples")
    
    hist(transect_data_natura$Effort, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Transect length",
         main = "Histogram of transect lengths")
    
    hist(transect_data_natura$PatchDensity, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Patch density",
         main = "Histogram of Natura patch density in transects")
    
    hist(transect_data_corine$PatchDensity, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Patch density",
         main = "Histogram of Corine patch density in transects")
    
    hist(transect_data_natura$SimpsonsDiversity, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Simpson's diversity",
         main = "Histogram of Natura Simpson's diversity in transects")
    
    hist(transect_data_corine$SimpsonsDiversity, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Simpson's diversity",
         main = "Histogram of Corine Simpson's diversity in transects")
    
    hist(transect_data_natura$ShannonsDiversity, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Shannon's diversity",
         main = "Histogram of Natura Shannon's diversity in transects")
    
    hist(transect_data_corine$ShannonsDiversity, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Shannon's diversity",
         main = "Histogram of Corine Shannon's diversity in transects")
    
    hist(transect_data_natura$ScaledRichness, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Scaled richness",
         main = "Histogram of Natura scaled richness in transects")
    
    hist(transect_data_corine$ScaledRichness, 
         probability = TRUE, 
         col = "lightgray",
         xlab = "Scaled richness",
         main = "Histogram of Corine scaled richness in transects")
    
    
    # Spatial plots for variables
    
    # Habitat type clusters spatially
    natura_clusters <- get_clusters(fractions_natura, 5)
    plot(coordinates[names(natura_clusters), ]$x,
         coordinates[names(natura_clusters), ]$y, 
         cex =  1, 
         col = natura_clusters,
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Natura clusters spatial pattern")
    corine_clusters <- get_clusters(fractions_corine, 5)
    plot(coordinates[names(corine_clusters), ]$x,
         coordinates[names(corine_clusters), ]$y, 
         cex =  1, 
         col = corine_clusters,
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Corine clusters spatial pattern")
    # Cluster value correlation
    cluster_correlation <- cor(natura_clusters, corine_clusters)
    plot(coordinates[names(corine_clusters), ]$x,
         coordinates[names(corine_clusters), ]$y, 
         cex =  1, 
         col = natura_clusters,
         pch = corine_clusters,
         xlab = "X", 
         ylab = "Y", 
         main = sprintf("Both (color: natura, shape: corine), correlation = %s", 
                        round(cluster_correlation, 2)))
    
    
    # Diversities spatially
    plot(coordinates[transect_data_natura$Transect, ]$x,
         coordinates[transect_data_natura$Transect, ]$y, 
         cex =  transect_data_natura$PatchDensity * 0.003, 
         col = "Black",
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Natura patch density")
    plot(coordinates[transect_data_corine$Transect, ]$x,
         coordinates[transect_data_corine$Transect, ]$y, 
         cex =  transect_data_corine$PatchDensity * 0.003, 
         col = "Black",
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Corine patch density")
    plot(coordinates[transect_data_natura$Transect, ]$x,
         coordinates[transect_data_natura$Transect, ]$y, 
         cex =  transect_data_natura$SimpsonsDiversity, 
         col = "Black",
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Natura Simpson's diversity")
    plot(coordinates[transect_data_corine$Transect, ]$x,
         coordinates[transect_data_corine$Transect, ]$y, 
         cex =  transect_data_corine$SimpsonsDiversity, 
         col = "Black",
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Corine Simpson's diversity")
    plot(coordinates[transect_data_natura$Transect, ]$x,
         coordinates[transect_data_natura$Transect, ]$y, 
         cex =  transect_data_natura$ShannonsDiversity, 
         col = "Black",
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Natura Shannon's diversity")
    plot(coordinates[transect_data_corine$Transect, ]$x,
         coordinates[transect_data_corine$Transect, ]$y, 
         cex =  transect_data_corine$ShannonsDiversity, 
         col = "Black",
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Corine Shannon's diversity")
    plot(coordinates[transect_data_natura$Transect, ]$x,
         coordinates[transect_data_natura$Transect, ]$y, 
         cex =  transect_data_natura$ScaledRichness * 0.05, 
         col = "Black",
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Natura scaled richness")
    plot(coordinates[transect_data_corine$Transect, ]$x,
         coordinates[transect_data_corine$Transect, ]$y, 
         cex =  transect_data_corine$ScaledRichness * 0.05, 
         col = "Black",
         pch = 21,
         xlab = "X", 
         ylab = "Y", 
         main = "Corine scaled richness")
    
    
    
    # Diversity correlations
    correlation_columns <- c("PatchDensity", "SimpsonsDiversity", "ShannonsDiversity", "ScaledRichness", "Cluster")
    within_natura <- cor(transect_data_natura[, correlation_columns], use = "complete.obs")
    within_corine <- cor(transect_data_corine[, correlation_columns], use = "complete.obs")
    between_datasets <- matrix(NA, 
                               nrow = length(correlation_columns), 
                               ncol = length(correlation_columns), 
                               dimnames = list(correlation_columns, correlation_columns))
    for (var1 in correlation_columns) {
        for (var2 in correlation_columns) {
            between_datasets[var1, var2] <- cor(transect_data_natura[,var1], 
                                                transect_data_corine[,var2], 
                                                use = "complete.obs")
        }
    }
    diversity_correlations <- list(
        "Diversity correlation within natura" = within_natura,
        "Diversity correlation within corine" = within_corine,
        "Diversity correlation between datasets" = between_datasets
    )
    for (name in names(diversity_correlations)) {
        pheatmap(diversity_correlations[[name]],
                 main = name,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 display_numbers = TRUE)
    }
    
    # Habitat type correlations
    pheatmap(cor(fractions_natura),
             main = "Natura type correlations",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             display_numbers = TRUE)
    
    pheatmap(cor(fractions_corine),
             main = "Corine type correlations",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             fontsize = 4,
             display_numbers = TRUE)
    
    
    
    # Temperatures spatially
    colors <- colorRampPalette(c("blue", "white", "red"))
    number_of_colors <- 100 
    scaled_temperatures <- scale_to_range(transect_data_natura$Temperature, 
                                          1, 
                                          number_of_colors)
    temperature_colors <- colors(number_of_colors)[round(scaled_temperatures)]
    plot(coordinates[transect_data_natura$Transect, ]$x,
         coordinates[transect_data_natura$Transect, ]$y, 
         cex =  1, 
         col = temperature_colors,
         pch = 20,
         xlab = "X", 
         ylab = "Y", 
         main = "Transect mean spring temperatures over all years")
    text(x = coordinates[transect_data_natura$Transect, ]$x,
         y = coordinates[transect_data_natura$Transect, ]$y,
         labels = round(transect_data_natura$Temperature, 1),  
         pos = 3, 
         col = "black", 
         cex = 0.5)
    
    
    # Yearly mean temperatures
    yearly_mean_temperatures <- aggregate(Temperature ~ Year, 
                                          data = data.frame(Temperature = natura$Temperature,
                                                            Year = spatiotemporal_context$Year),
                                          FUN = mean)
    barplot(
        height = yearly_mean_temperatures$Temperature,  
        names.arg = yearly_mean_temperatures$Year,     
        xlab = "Year",
        ylab = "Temperature (mean of April and May)",
        main = "Yearly mean spring temperatures over all transects")
    
    

    dev.off()
    
}




# SCRIPT STARTS 


# LOAD RAW DATA
load(file = file.path(dir_data, "occurrence_raw.RData")) 
load(file = file.path(dir_data, "abundance_raw.RData")) 
load(file = file.path(dir_data, "spatiotemporal_context_raw.RData"))
load(file = file.path(dir_data, "transect_coordinates.RData"))
load(file = file.path(dir_data, "env_data_natura_raw.RData")) 
load(file = file.path(dir_data, "env_data_corine_raw.RData"))
load(file = file.path(dir_data, "fractions_natura.RData")) 
load(file = file.path(dir_data, "fractions_corine.RData"))


explore_bird_data(occurrence, abundance, spatiotemporal_context, transect_coordinates, "raw")
explore_habitat_data(env_data_natura, 
                     env_data_corine, 
                     fractions_natura,
                     fractions_corine,
                     spatiotemporal_context, 
                     transect_coordinates, 
                     "raw")
