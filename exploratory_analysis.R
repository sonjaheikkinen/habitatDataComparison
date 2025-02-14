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
                                                                   iter.max = 15)
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

plot_silhouette_scores <- function(data, data_name) {
    max_clusters <- 15
    dist_data <- dist(data)
    clusters <- hclust(dist_data, method = "complete")
    silhouette_scores <- sapply(2:max_clusters,
                               function(cluster_amount) {
                                   clusters <- cutree(clusters, cluster_amount)
                                   silhouette_scores <- silhouette(clusters, dist_data)
                                   mean_silhouette_score <- mean(silhouette_scores[,3])
                                   return(mean_silhouette_score)
                               })
    plot(2:max_clusters,
         silhouette_scores,
         frame = FALSE,
         type = "b",
         xlab = "Number of clusters",
         ylab = "Average silhouette scores",
         main = sprintf("Silhouette scores for %s", data_name))
}



scale_to_range <- function(values, new_min, new_max) {
    old_min <- min(values)
    old_max <- max(values)
    scaled_values <- new_min + (values - old_min) * (new_max - new_min) / (old_max - old_min)
    return(scaled_values)
}



plot_spatially <- function(values, 
                           sizes, 
                           coord_indices, 
                           coordinates, 
                           title, 
                           colors = "Black",
                           pch = 21) {
    plot(coordinates[coord_indices, ]$x,
         coordinates[coord_indices, ]$y, 
         cex = sizes, 
         pch = pch,
         col = colors,
         xlab = sprintf("X \n[%s, %s], mean %s, median %s, sd %s",
                               round(min(values), 2), 
                               round(max(values), 2),
                               round(mean(values), 2), 
                               round(median(values), 2),
                               round(sd(values), 2)), 
         ylab = "Y", 
         main = title)
}

plot_barplot <- function(values, names, title) {
    barplot(values, 
            names.arg = names,
            horiz = TRUE,
            las = 2,
            xlab = sprintf("[%s, %s], mean %s, median %s, sd %s",
                           round(min(values), 2), 
                           round(max(values), 2),
                           round(mean(values), 2), 
                           round(median(values), 2),
                           round(sd(values), 2)),
            main = title)
}



# EXPLORATORY ANALYSIS FUNCTIONS

explore_bird_data <- function(occurrence, abundance, spatiotemporal_context, coordinates, type) {
    
    pdf(file.path(dir_results, sprintf("exploratory_analysis_observations_%s.pdf", type)))
    
    
    sample_total_abundance <- rowSums(abundance, na.rm = TRUE)
    sample_species_richness <- rowSums(occurrence, na.rm = TRUE)
    sample_simpson_diversity <- diversity(abundance, index = "simpson")
    sample_shannon_diversity <- diversity(abundance, index = "shannon")
    
    
    # SPATIAL AND TEMPORAL PLOTS FOR TOTAL ABUNDANCE, SPECIES RICHNESS, DIVERSITIES
    
    
    # Spatial overview
    transect_mean_total_abundance <- aggregate(Abundance ~ Transect,
                                               data = data.frame(Abundance = sample_total_abundance,
                                                                 spatiotemporal_context),
                                               FUN = mean)
    transect_mean_species_richness <- aggregate(Richness ~ Transect,
                                               data = data.frame(Richness = sample_species_richness,
                                                                 spatiotemporal_context),
                                               FUN = mean)
    transect_mean_simpson_diversity <- aggregate(Diversity ~ Transect,
                                               data = data.frame(Diversity = sample_simpson_diversity,
                                                                 spatiotemporal_context),
                                               FUN = mean)
    transect_mean_shannon_diversity <- aggregate(Diversity ~ Transect,
                                                 data = data.frame(Diversity = sample_shannon_diversity,
                                                                   spatiotemporal_context),
                                                 FUN = mean)
    plot_spatially(transect_mean_total_abundance$Abundance,
                   transect_mean_total_abundance$Abundance * 0.005,
                   transect_mean_total_abundance$Transect, 
                   coordinates,
                   "Mean of sample total abundances for each transect")
    
    plot_spatially(transect_mean_species_richness$Richness,
                   transect_mean_species_richness$Richness * 0.05,
                   transect_mean_species_richness$Transect,
                   coordinates,
                   "Mean of sample species richnesses for each transect")
    scaled_simpson_diversities <- scale(transect_mean_simpson_diversity$Diversity)
    plot_spatially(transect_mean_simpson_diversity$Diversity,
                   (scaled_simpson_diversities + abs(min(scaled_simpson_diversities))) * 0.3,
                   transect_mean_simpson_diversity$Transect,
                   coordinates,
                   "Mean of sample Simpson's diversities for each transect")
    scaled_shannon_diversities <- scale(transect_mean_shannon_diversity$Diversity)
    plot_spatially(transect_mean_shannon_diversity$Diversity,
                   (scaled_shannon_diversities + abs(min(scaled_shannon_diversities))) * 0.3,
                   transect_mean_shannon_diversity$Transect,
                   coordinates,
                   "Mean of sample Shannon's diversities for each transect")
    
    
    # Temporal overview
    year_mean_total_abundance <- aggregate(Abundance ~ Year,
                                           data = data.frame(Abundance = sample_total_abundance,
                                                             spatiotemporal_context),
                                           FUN = mean)
    year_mean_species_richness <- aggregate(Richness ~ Year,
                                            data = data.frame(Richness = sample_species_richness,
                                                              spatiotemporal_context),
                                            FUN = mean)
    year_mean_simpson_diversity <- aggregate(Diversity ~ Year,
                                             data = data.frame(Diversity = sample_simpson_diversity,
                                                               spatiotemporal_context),
                                             FUN = mean)
    year_mean_shannon_diversity <- aggregate(Diversity ~ Year,
                                             data = data.frame(Diversity = sample_shannon_diversity,
                                                               spatiotemporal_context),
                                             FUN = mean)
    plot_barplot(year_mean_total_abundance$Abundance, 
                 year_mean_total_abundance$Year,
                 "Mean of sample total abundances for each year")
    plot_barplot(year_mean_species_richness$Richness, 
                 year_mean_species_richness$Year,
                 "Mean of sample species richness for each year")
    plot_barplot(year_mean_simpson_diversity$Diversity, 
                 year_mean_simpson_diversity$Year,
                 "Mean of sample Simpson's diversities for each year")
    plot_barplot(year_mean_shannon_diversity$Diversity, 
                 year_mean_shannon_diversity$Year,
                 "Mean of sample Shannon's diversities for each year")
    
    
    # Yearly plots
    
    
    
    
    
    
    # Total abundance in samples
    hist(sample_total_abundance,
         main = "Histogram of total abundance in samples",
         xlab = sprintf("Abundance, min %s, max %s", min(sample_total_abundance), max(sample_total_abundance)),
         ylab = sprintf("Frequency, %s total samples", nrow(abundance)))
    
    # Total species richness in samples
    hist(sample_species_richness,
         main = "Histogram of species richness in samples",
         xlab = sprintf("Occurrence, min %s, max %s", min(sample_species_richness), max(sample_species_richness)),
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
    sample_total_abundance <- matrix(0,
                                      nrow = length(unique(spatiotemporal_context$Year)),
                                      ncol = length(unique(spatiotemporal_context$Transect)))
    rownames(sample_total_abundance) <- rev(sort(unique(spatiotemporal_context$Year)))
    colnames(sample_total_abundance) <- unique(spatiotemporal_context$Transect)
    for (row in 1:nrow(spatiotemporal_context)) {
        year <- spatiotemporal_context[row, ]$Year
        transect <- spatiotemporal_context[row, ]$Transect
        sample <- sprintf("%s%s", year, transect)
        year_index <- match(year, rownames(sample_total_abundance))
        transect_index <- match(transect, colnames(sample_total_abundance))
        row_abundances <- abundance[row, ]
        row_abundances[is.na(row_abundances)] <- 0 
        sample_total_abundance[year_index, transect_index] <- sum(row_abundances)
    }
    pheatmap(sample_total_abundance,
             cluster_rows = FALSE,
             main = "Total abundance in each sample")
    
    
    # Spatial variation in transect total abundance 
    # Sample total abundances: sum of species abundances in each sample
    # Transect mean total abundances: mean of sample total abundances from that transect
    transect_mean_total_abundances <- apply(sample_total_abundance, 
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
    year_mean_total_abundances <- apply(sample_total_abundance, 
                                   1, 
                                   function(row) {
                                       non_zero_values <- row[row != 0]
                                       if (length(non_zero_values) > 0) {
                                           mean(non_zero_values)
                                       } else {
                                           NA  # Return NA if there are no non-zero values
                                       }
                                   })
    names(year_mean_total_abundances) <- rownames(sample_total_abundance)
    barplot(year_mean_total_abundances,
            horiz = TRUE,
            las = 1,
            main = "Year mean total abundances")
    
    
    
    # Spatial variation in species composition
    
    # Create spatial occurrence and abundance data
    transect_occurrence_data_list <- list()
    transect_abundance_data_list <- list()
    for (transect in unique(spatiotemporal_context$Transect)) {
        transect_indices <- which(spatiotemporal_context$Transect == transect)
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
    transect_mean_simpson_diversities <- c()
    transect_mean_shannon_diversities <- c()
    for (transect in unique(spatiotemporal_context)$Transect) {
        transect_indices <- which(spatiotemporal_context$Transect == transect)
        transect_mean_simpson_diversities[transect] <- mean(sample_simpson_diversity[transect_indices])
        transect_mean_shannon_diversities[transect] <- mean(sample_shannon_diversity[transect_indices])
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
    #plot_cluster_elbow(transect_occurrence_data, "transect species occurrences")
    plot_silhouette_scores(transect_occurrence_data, "transect species occurrences")
    
    # Are there truly different transects based on transect abundance
    #plot_cluster_elbow(transect_abundance_data, "transect species abundances")
    plot_silhouette_scores(transect_abundance_data, "transect species abundances")
    
    
    dev.off()
    
}


# FUNCTIONS | HABITAT DATA

plot_distribution <- function(data, xlab, main) {
    hist(data, 
         xlab = sprintf("%s, \n[%s, %s], mean %s, median %s, sd %s",
                        xlab, 
                        round(min(data), 2), 
                        round(max(data), 2),
                        round(mean(data), 2), 
                        round(median(data), 2),
                        round(sd(data), 2)),
         main = main)
}


create_pca_plots <- function(pca_results, title) {
    # For each pca, calculate the amount of variation it accounts for (eigenvalue)
    eigenvalues <- pca_results$sdev^2
    # Transform the actual values for percentages for easier understanding
    percentage_of_variation_accounted_for_by_pca_axes <- round(eigenvalues / sum(eigenvalues) * 100, 1)
    plot_pca_scree(percentage_of_variation_accounted_for_by_pca_axes, title)
    plot_pca_eigenvalues(eigenvalues, title)
    plot_pca_clusters(pca_results, 
                      percentage_of_variation_accounted_for_by_pca_axes, 
                      title)
    plot_pca_loadings(pca_results, title, pca_number = 1)
    plot_pca_loadings(pca_results, title, pca_number = 2)
}

plot_pca_scree <- function(percentage_of_variation_accounted_for_by_pca_axes,
                           title) {
    
    cumulative_variation <- cumsum(percentage_of_variation_accounted_for_by_pca_axes)
    
    plot(1:length(percentage_of_variation_accounted_for_by_pca_axes), 
         percentage_of_variation_accounted_for_by_pca_axes, 
         type = "b", 
         col = "darkgreen", 
         pch = 19, 
         ylim = c(0, 100), 
         xlab = "Principal Component", 
         ylab = "Percentage of Variation", 
         main = sprintf("%s Scree Plot", title))
    
    lines(1:length(cumulative_variation), 
          cumulative_variation, 
          type = "b", 
          col = "blue", 
          pch = 19)
    
    abline(h = seq(0, 100, by = 10), col = "gray", lty = "dotted")
    
    legend("topright", 
           legend = c("Individual Variation", "Cumulative Variation"), 
           col = c("darkgreen", "blue"), 
           pch = 19, 
           lty = 1)
    
}

plot_pca_eigenvalues <- function(eigenvalues, title) {
    plot(1:length(eigenvalues), 
         eigenvalues, 
         type = "b", 
         col = "blue", 
         pch = 19, 
        xlab = "Principal Component", 
        ylab = "Eigenvalue", 
        main = sprintf("%s Eigenvalues", title))
    
    abline(h = 1, col = "black", lty = "dashed")
    
    legend("topright", 
           legend = c("Eigenvalues", "Threshold (y=1)"), 
           col = c("blue", "black"), 
           pch = c(19, NA), 
           lty = c(1, 2)
    )
}


plot_pca_clusters <- function(pca_results, pca_variation_percentage, title) {
    pca_plot_data <- data.frame(Transect = rownames(pca_results$x), # Transect numbers
                                X = pca_results$x[,1], # pc1 coordinate of transect
                                Y = pca_results$x[,2]) # pc2 coordinate of transect
    print(ggplot(data = pca_plot_data, aes(x = X, y = Y, label = Transect)) +
              geom_text(size = 2) +
              xlab(paste("PC1 - ", pca_variation_percentage[1], " %", sep = "")) +
              ylab(paste("PC2 - ", pca_variation_percentage[2], " %", sep = "")) +
              theme_bw() +
              ggtitle(sprintf("%s PCA", title)))
}


plot_pca_loadings <- function(pca_results, title, pca_number) {
    loading_scores <- pca_results$rotation[,pca_number]
    barplot(sort(loading_scores), 
            horiz = TRUE, 
            las = 2,
            main = sprintf("%s, loadings for PC1", title))
}



explore_habitat_data <- function(natura, 
                                 corine, 
                                 fractions_natura,
                                 fractions_corine,
                                 pca_results_natura,
                                 pca_results_corine,
                                 spatiotemporal_context, 
                                 coordinates, 
                                 occurrence,
                                 abundance,
                                 type) {
    
    
    pdf(file.path(dir_results, sprintf("exploratory_analysis_habitats_%s.pdf", type)))
    
    
    create_pca_plots(pca_results_natura, "Natura")
    create_pca_plots(pca_results_corine, "Corine")
    
    
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
    

    # Plot histograms of natura and corine types in transects
    # Also plot the histogram of just those values that aren't zeros
    for (column in colnames(fractions_natura)) {
        par(mfrow = c(2, 1))
        column_data <- fractions_natura[, column]
        column_data_without_zeros <- column_data[column_data != 0]
        non_zeros <- sum(column_data != 0)
        plot_distribution(fractions_natura[,column], 
                          xlab = sprintf("%s", column),
                          main = sprintf("Natura: Histogram of %s in transects, \nnon-zeros: %s/%s = %s", 
                                         column, 
                                         non_zeros, 
                                         length(column_data), 
                                         round(non_zeros / length(column_data), 2)))
        if (length(column_data_without_zeros) > 0) {
            plot_distribution(column_data_without_zeros, 
                              xlab = sprintf("%s", column),
                              main = sprintf("Natura: Histogram of %s without zeros", column))
        }
    }
    for (column in colnames(fractions_corine)) {
        par(mfrow = c(2, 1))
        column_data <- fractions_corine[, column]
        column_data_without_zeros <- column_data[column_data != 0]
        non_zeros <- sum(column_data != 0)
        plot_distribution(fractions_corine[,column],
                          xlab = sprintf("%s", column),
                          main = sprintf("Corine: Histogram of %s in transects, \nnon-zeros: %s/%s = %s",
                                         column, 
                                         non_zeros, 
                                         length(column_data), 
                                         round(non_zeros / length(column_data), 2)))
        if (length(column_data_without_zeros) > 0) {
            plot_distribution(column_data_without_zeros,
                              xlab = sprintf("%s", column),
                              main = sprintf("Corine: Histogram of %s without zeros", column))
        }
    }
    
    par(mfrow=c(1,1))
    
    plot_distribution(natura$Temperature,
                      xlab = "Temperature (april + may)",
                      main = "Histogram of transect spring temperatures in samples")
    
    plot_distribution(natura$Effort,
                      xlab = "Transect length",
                      main = "Histogram of transect lengths in samples")
    
    variables_to_plot <- c("PatchDensity", 
                           "SimpsonsDiversity",
                           "ShannonsDiversity",
                           "ScaledRichness")
    
    for (column in variables_to_plot) {
        par(mfrow = c(2, 1))
        plot_distribution(transect_data_natura[,column],
                          xlab = column,
                          main = sprintf("Histogram of Natura %s in transects",
                                         column))
        plot_distribution(transect_data_corine[,column],
                          xlab = column,
                          main = sprintf("Histogram of Corine %s in transects",
                                         column))
    }
    
    par(mfrow=c(1,1))
    
    
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
    
    
    
    # Correlations
    correlation_columns <- c("PatchDensity", 
                             "SimpsonsDiversity", 
                             "ShannonsDiversity", 
                             "ScaledRichness", 
                             "Cluster",
                             "Temperature")
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
        "Correlation within natura" = within_natura,
        "Correlation within corine" = within_corine,
        "Correlation between datasets" = between_datasets
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
    
    
    
    # Temperatures spatially (means per transect)
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
    
    # Temperatures spatially (per year)
    for (year in unique(spatiotemporal_context$Year)) {
        year_indices <- which(spatiotemporal_context$Year == year)
        data_for_year <- data.frame(spatiotemporal_context[year_indices, ],
                                         natura[year_indices, ])
        colors <- colorRampPalette(c("blue", "white", "red"))
        number_of_colors <- 100 
        scaled_temperatures <- scale_to_range(data_for_year$Temperature, 
                                              1, 
                                              number_of_colors)
        temperature_colors <- colors(number_of_colors)[round(scaled_temperatures)]
        plot(coordinates[data_for_year$Transect, ]$x,
             coordinates[data_for_year$Transect, ]$y, 
             cex =  1, 
             col = temperature_colors,
             pch = 20,
             xlab = "X", 
             ylab = "Y", 
             main = sprintf("Transect mean spring temperatures for %s", year))
        text(x = coordinates[data_for_year$Transect, ]$x,
             y = coordinates[data_for_year$Transect, ]$y,
             labels = round(data_for_year$Temperature, 1),  
             pos = 3, 
             col = "black", 
             cex = 0.5)
        
    }
    
    
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
    
    
    # Elbow curvers and silhouettes
    
    #plot_cluster_elbow(fractions_natura, "Natura fractions")
    plot_silhouette_scores(fractions_natura, "Natura fractions")
    #plot_cluster_elbow(fractions_corine, "Corine fractions")
    plot_silhouette_scores(fractions_corine, "Corine fractions")
    
    diversity_columns <- c("PatchDensity", 
                           "SimpsonsDiversity", 
                           "ScaledRichness")
    #plot_cluster_elbow(transect_data_natura[,diversity_columns], "Natura diversity")
    plot_silhouette_scores(transect_data_natura[,diversity_columns], "Natura diversity")
    #plot_cluster_elbow(transect_data_corine[,diversity_columns], "Corine diversity")
    plot_silhouette_scores(transect_data_corine[,diversity_columns], "Corine diversity")
    
    #plot_cluster_elbow(transect_data_natura[,c("Temperature")], "Natura temperature")
    plot_silhouette_scores(transect_data_natura[,c("Temperature")], "Natura temperature")
    #plot_cluster_elbow(transect_data_corine[,c("Temperature")], "Corine temperature")
    plot_silhouette_scores(transect_data_corine[,c("Temperature")], "Corine temperature")
    
    
    combination_columns_natura <- c(colnames(fractions_natura), diversity_columns, "Temperature")
    combination_columns_corine <- c(colnames(fractions_corine), diversity_columns, "Temperature")
    #plot_cluster_elbow(transect_data_natura[,combination_columns_natura], "Natura combined")
    plot_silhouette_scores(transect_data_natura[,combination_columns_natura], "Natura combined")
    #plot_cluster_elbow(transect_data_corine[,combination_columns_corine], "Corine combined")
    plot_silhouette_scores(transect_data_corine[,combination_columns_corine], "Corine combined")
    
    
    
    # Plot species abundances on environmental gradients
    abundance_clusters <- get_clusters(as.data.frame(abundance), 5)
    
    # For all of these: spatial plot, temporal plot, linear model?
    
    # Create fractions datasets
    # Natura
    natura_cluster_abundances <- list()
    natura_cluster_richnesses <- list()
    natura_cluster_diversity <- list()
    natura_cluster_species_clusters <- list()
    for (cluster in unique(natura$Cluster)) {
        cluster_indices <- which(natura$Cluster == cluster)
        cluster_abundance_rows <- abundance[cluster_indices, ]
        cluster_abundances <- rowSums(cluster_abundance_rows)
        cluster_richnesses <- apply(cluster_abundance_rows,
                                    1,
                                    function(row) {
                                        return(sum(row > 0))
                                    })
        cluster_diversities <- diversity(cluster_abundance_rows, index = "simpson")
        cluster_species_clusters <- abundance_clusters[cluster_indices]
        natura_cluster_abundances[cluster] <- list(cluster_abundances)
        natura_cluster_richnesses[cluster] <- list(cluster_richnesses)
        natura_cluster_diversity[cluster] <- list(cluster_diversities)
        natura_cluster_species_clusters[cluster] <- list(cluster_species_clusters)
    }
    # Corine
    corine_cluster_abundances <- list()
    corine_cluster_richnesses <- list()
    corine_cluster_diversity <- list()
    corine_cluster_species_clusters <- list()
    for (cluster in unique(corine$Cluster)) {
        cluster_indices <- which(corine$Cluster == cluster)
        cluster_abundance_rows <- abundance[cluster_indices, ]
        cluster_abundances <- rowSums(cluster_abundance_rows)
        cluster_richnesses <- apply(cluster_abundance_rows,
                                    1,
                                    function(row) {
                                        return(sum(row > 0))
                                    })
        cluster_diversities <- diversity(cluster_abundance_rows, index = "simpson")
        cluster_species_clusters <- abundance_clusters[cluster_indices]
        corine_cluster_abundances[cluster] <- list(cluster_abundances)
        corine_cluster_richnesses[cluster] <- list(cluster_richnesses)
        corine_cluster_diversity[cluster] <- list(cluster_diversities)
        corine_cluster_species_clusters[cluster] <- list(cluster_species_clusters)
    }
     
    # Also glms for all these?
    
    # Fractions vs abundance
    boxplot(natura_cluster_abundances,
            xlab = "Natura cluster",
            ylab = "Abundance",
            main = "Sample total abundance distributions for each Natura Cluster")
    boxplot(corine_cluster_abundances,
            xlab = "Corine cluster",
            ylab = "Abundance",
            main = "Sample total abundance distributions for each Corine Cluster")
    
    # Fractions vs richness
    boxplot(natura_cluster_richnesses,
            xlab = "Natura cluster",
            ylab = "Species richnesses",
            main = "Sample species richnesses for each Natura Cluster")
    boxplot(corine_cluster_richnesses,
            xlab = "Corine cluster",
            ylab = "Species richnesses",
            main = "Sample species richnesses for each Corine Cluster")
    
    # Fractions vs diversity
    boxplot(natura_cluster_diversity,
            xlab = "Natura cluster",
            ylab = "Simpson's diversity",
            main = "Sample Simpson's diversities for each Natura Cluster")
    boxplot(corine_cluster_diversity,
            xlab = "Corine cluster",
            ylab = "Simpson's diversity",
            main = "Sample Simpson's diversities for each Corine Cluster")
    
    # Fractions vs species cluster
    # Should this be boxplot?
    
    # Diversity vs abundance
    
    # Diversity vs richness
    
    # Diversity vs diversity
    
    # Diversity vs. species cluster

    # Glms for all these also?
    
    # Visualize what natura and corine clusters mean in practice (habitat fractions)    
    
    
    
    dev.off()
    
}


explore_phylogeny_data <- function(tree, abundance, type) {
    
    pdf(file.path(dir_results, sprintf("exploratory_analysis_phylogeny_%s.pdf", type)))
    
    # First group the species to order level
    order_names <- c("Sorsalinnut",
                     "Kanalinnut",
                     "Kuikkalinnut",
                     "Rantalinnut",
                     "Kiitäjälinnut",
                     "Pöllölinnut",
                     "Tikkalinnut",
                     "Päiväpetolinnut",
                     "Varpuslinnut",
                     "Jalohaukkalinnut",
                     "Käkilinnut")
    groups <- cutree(tree, 12) # Jaa lahkotasolle (manuaalisesti selvitetty)
    
    # Then prune tree to those in abundance
    species_list <- colnames(abundance)
    species_to_remove <- setdiff(tree$tip.label, species_list)
    tree <- drop.tip(tree, species_to_remove)
    groups <- groups[setdiff(names(groups), species_to_remove)]
    
    # Plot phylogenetic tree
    plot(tree,
         cex = 0.3)
    
    # Plot phylogenetic correlation matrix
    correlation_matrix <- vcv(tree)
    pheatmap(correlation_matrix,
             fontsize_row = 4,
             fontsize_col = 4)
    

    group_total_abundances <- sapply(unique(groups), 
                                     function(group) {
                                         group_indices <- which(groups == group)
                                         group_species <- names(groups)[group_indices]
                                         print(group)
                                         print(group_species)
                                         group_abundance <- sum(abundance[,group_species])
                                         return(group_abundance)
                                     })
    names(group_total_abundances) <- order_names
    group_species_richnesses <- sapply(unique(groups), 
                                       function(group) {
                                           group_indices <- which(groups == group)
                                           return(length(group_indices))
                                       })
    names(group_species_richnesses) <- order_names
    
    # Adjust margin to fit species names to the left
    old_par <- par(no.readonly = TRUE) 
    par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))
    
    barplot(group_species_richnesses,
            horiz = TRUE,
            las = 1,
            main = "Number of species per order")
    barplot(group_total_abundances,
            horiz = TRUE,
            las = 1,
            main = "Total abundance per order (lahko)")
    barplot(group_total_abundances / group_species_richnesses,
            horiz = TRUE,
            las = 1,
            main = "Mean abundance per order (total / number of species)")
  
    par(old_par)
    
    
    dev.off()
    
}


# FUNCTIONS | TRAIT DATA

count_trait_total_abundance <- function(trait_value, 
                                        column, 
                                        trait_data, 
                                        species_total_abundances) {
    trait_value_indices <- which(trait_data[,column] == trait_value)
    species_with_trait_value <- trait_data[trait_value_indices, ]$Species
    trait_value_abundance <- sum(species_total_abundances[species_with_trait_value])
    return(trait_value_abundance)
}

explore_trait_data <- function(traits, abundance, type) {
    
    pdf(file.path(dir_results, sprintf("exploratory_analysis_traits_%s.pdf", type)))
    
    barplot(table(trait_data$Feeding),
            main = "Feeding types of species",
            ylab = "Frequency")
    
    barplot(table(trait_data$Mig),
            main = "Migration types of species",
            ylab = "Frequency")
    
    barplot(table(trait_data$Habitat),
            main = "Habitat types of species",
            ylab = "Frequency")
    
    hist(trait_data$Mass,
         main = "Histogram of mass of species")
    
    # Traits weighted by abundances in the data
    species_total_abundances <- colSums(abundance)
    feeding_types_in_data <- sapply(unique(trait_data$Feeding),
                                    function(trait_value) {
                                        return(count_trait_total_abundance(trait_value,
                                                                           "Feeding",
                                                                           trait_data, 
                                                                           species_total_abundances))
                                    })
    migration_types_in_data <- sapply(unique(trait_data$Mig),
                                    function(trait_value) {
                                        return(count_trait_total_abundance(trait_value,
                                                                           "Mig",
                                                                           trait_data, 
                                                                           species_total_abundances))
                                    })
    habitat_types_in_data <- sapply(unique(trait_data$Habitat),
                                    function(trait_value) {
                                        return(count_trait_total_abundance(trait_value,
                                                                           "Habitat",
                                                                           trait_data, 
                                                                           species_total_abundances))
                                    })
    masses_in_data <- unlist(sapply(trait_data$Species,
                             function(species) {
                                 species_mass <- trait_data[trait_data$Species == species,]$Mass
                                 species_total_abundance <- species_total_abundances[species]
                                 return(rep(species_mass, species_total_abundance))
                             }))
    names(masses_in_data) <- NULL
    barplot(feeding_types_in_data,
         main = "Total abundance for each feeding type in data")
    barplot(migration_types_in_data,
            main = "Total abundance for each migration type in data")
    barplot(habitat_types_in_data,
            main = "Total abundance for each habitat type in data")
    hist(masses_in_data,
         main = "Total abundance of all masses")
    hist(masses_in_data[masses_in_data >= 1000],
         main = "Histogram of masses over 1000 g")
    hist(masses_in_data[masses_in_data < 1000 & masses_in_data >= 100],
         main = "Histogram of masses under 1000 g but over 100 g")
    hist(masses_in_data[masses_in_data < 100],
         main = "Histogram of masses under 100 g")
    

    
    
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
load(file = file.path(dir_data, "fractions_natura_raw.RData")) 
load(file = file.path(dir_data, "pca_results_natura_raw.RData"))
load(file = file.path(dir_data, "pca_results_corine_raw.RData"))
load(file = file.path(dir_data, "fractions_corine_raw.RData"))
load(file = file.path(dir_data, "phylogeny_data_raw.RData"))
load(file = file.path(dir_data, "trait_data_raw.RData"))


explore_bird_data(occurrence, 
                  abundance, 
                  spatiotemporal_context, 
                  transect_coordinates, 
                  "raw")
explore_habitat_data(env_data_natura, 
                     env_data_corine, 
                     fractions_natura,
                     fractions_corine,
                     pca_results_natura,
                     pca_results_corine,
                     spatiotemporal_context, 
                     transect_coordinates, 
                     occurrence,
                     abundance,
                     "raw")
explore_phylogeny_data(phylogeny_data, 
                       abundance,
                       "raw")
explore_trait_data(trait_data,
                   abundance, 
                   "raw")



# LOAD SELECTED DATA
load(file = file.path(dir_data, "occurrence.RData")) 
load(file = file.path(dir_data, "abundance.RData")) 
load(file = file.path(dir_data, "spatiotemporal_context.RData"))
load(file = file.path(dir_data, "transect_coordinates.RData"))
load(file = file.path(dir_data, "env_data_natura.RData")) 
load(file = file.path(dir_data, "env_data_corine.RData"))
load(file = file.path(dir_data, "fractions_natura.RData")) 
load(file = file.path(dir_data, "fractions_corine.RData"))
load(file = file.path(dir_data, "pca_results_natura.RData"))
load(file = file.path(dir_data, "pca_results_corine.RData"))
load(file = file.path(dir_data, "phylogeny_data.RData"))
load(file = file.path(dir_data, "trait_data.RData"))


explore_bird_data(occurrence, 
                  abundance, 
                  spatiotemporal_context, 
                  transect_coordinates, 
                  "selected")
explore_habitat_data(env_data_natura, 
                     env_data_corine, 
                     fractions_natura,
                     fractions_corine,
                     pca_results_natura,
                     pca_results_corine,
                     spatiotemporal_context, 
                     transect_coordinates, 
                     occurrence,
                     abundance,
                     "selected")
explore_phylogeny_data(phylogeny_data, 
                       abundance,
                       "selected")
explore_trait_data(trait_data,
                   abundance, 
                   "selected")
