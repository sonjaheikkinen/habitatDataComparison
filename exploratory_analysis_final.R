# LOAD RAW DATA
observations <- read.csv(file.path(dir_data, "observations_not_year_limited.csv"),
                         sep = ";")

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
transects_shp <- vect(file.path(dir_data, "transects_selected.shp"))
natura_classification <- read_excel(file.path(dir_data, "Ylalappi_luokitus.xls"))
corine_classification <- read_excel(file.path(dir_data, "CorineMaanpeite2018Luokat.xls"))


# NAME MAPPING FOR HABITAT
names_natura <- data.frame(value = natura_classification$Value, name = natura_classification$NaturaTyyppi)
names_corine <- data.frame(value = corine_classification$Value, name = corine_classification$Level4Suo)

# TRANSECT NAMES FOR COORDINATE DATA
transect_coordinates$Transect <- rownames(transect_coordinates)




pdf(file.path(dir_results, sprintf("exploratory_analysis_final.pdf")))


# SPATIOTEMPORAL_CONTEXT OVERVIEW
years <- sort(unique(observations$Year))
transects <- unique(observations$SiteID)
samples <- matrix(0, nrow = length(transects), ncol = length(years))
colnames(samples) <- years
rownames(samples) <- transects
for (sample in 1:nrow(observations)) {
    year <- as.character(observations[sample, ]$Year)
    transect <- observations[sample, ]$SiteID
    samples[transect, year] <- 1
}
pheatmap(samples,
         color = c("lightblue", "navy"),
         cluster_cols = FALSE,
         legend_breaks = c(0, 1),
         main = "Transect visits")



# HABITAT RASTERS OVERVIEW
par(mfrow = c(1, 2))
plot(natura_raster,
     main = "Natura",
     legend = FALSE)
lines(transects_shp, col = "red")
corine_cropped <- crop(corine_raster, natura_raster)
plot(corine_cropped, 
     main = "Corine",
     legend = FALSE)
lines(transects_shp, col = "red")
title(main = "Natura and Corine types in study area, with included transects",
      outer = TRUE,
      line = -3)




# HABITAT TYPES OVERVIEW

# Frequency of values in the entire raster
natura_total_frequencies <- freq(natura_raster)
natura_total_frequencies$fraction <- natura_total_frequencies$count / sum(natura_total_frequencies$count)
natura_total_frequencies$region <- "Total"
corine_cropped <- crop(corine_raster, natura_raster)
corine_total_frequencies <- freq(corine_cropped)
corine_total_frequencies$fraction <- corine_total_frequencies$count / sum(corine_total_frequencies$count)
corine_total_frequencies$region <- "Total"

# Read lines and buffer
buffer_area <- buffer(transects_shp, width = buffer_width)

# Mask raster with buffer
natura_masked <- mask(natura_raster, buffer_area)
corine_masked <- mask(corine_cropped, buffer_area)

# Frequency within buffer
natura_buffer_frequencies <- freq(natura_masked)
natura_buffer_frequencies$fraction <- natura_buffer_frequencies$count / sum(natura_buffer_frequencies$count)
natura_buffer_frequencies$region <- "Buffer"
corine_buffer_frequencies <- freq(corine_masked)
corine_buffer_frequencies$fraction <- corine_buffer_frequencies$count / sum(corine_buffer_frequencies$count)
corine_buffer_frequencies$region <- "Buffer"

# Combine both data frames
natura_frequencies_combined <- rbind(natura_total_frequencies, natura_buffer_frequencies)
natura_frequencies_combined$value <- as.factor(natura_frequencies_combined$value)
corine_frequencies_combined <- rbind(corine_total_frequencies, corine_buffer_frequencies)
corine_frequencies_combined$value <- as.factor(corine_frequencies_combined$value)

# Add habitat type names
corine_frequencies_combined <- merge(
    corine_frequencies_combined, 
    names_corine[, c("value", "name")], 
    by.x = "value", 
    by.y = "value",
    all.x = TRUE
)

# Add habitat type names
natura_frequencies_combined <- merge(
    natura_frequencies_combined, 
    names_natura[, c("value", "name")], 
    by.x = "value", 
    by.y = "value",
    all.x = TRUE
)


natura_barplot <- ggplot(natura_frequencies_combined, aes(x = name, y = fraction, fill = region)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() + 
    labs(
        title = "1A Natura total fractions",
        x = "Habitat type",
        y = "Fraction",
        fill = "Region"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("Total" = "dodgerblue", "Buffer" = "darkblue"))

corine_barplot <- ggplot(corine_frequencies_combined, aes(x = name, y = fraction, fill = region)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(
        title = "1B Corine total fractions",
        x = "Land cover type",
        y = "Fraction",
        fill = "Region"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("Total" = "dodgerblue", "Buffer" = "darkblue"))





# Convert transect habitat fractions to long format
natura_fractions_long <- data.frame(
    variable = rep(colnames(fractions_natura), each = nrow(fractions_natura)),
    fraction = as.vector(as.matrix(fractions_natura))
)
corine_fractions_long <- data.frame(
    variable = rep(colnames(fractions_corine), each = nrow(fractions_corine)),
    fraction = as.vector(as.matrix(fractions_corine))
)


natura_boxplot <- ggplot(natura_fractions_long, aes(x = variable, y = fraction)) +
    geom_boxplot() +
    coord_flip() +
    labs(title = "2A Natura transect fractions", 
         x = "Habitat type", 
         y = "Fraction") +
    theme_minimal()

corine_boxplot <- ggplot(corine_fractions_long, aes(x = variable, y = fraction)) +
    geom_boxplot() +
    coord_flip() +
    labs(title = "2B Corine transect fractions", 
         x = "Land cover type", 
         y = "Fraction") +
    theme_minimal()


# Arrange all plots to grid
grid.arrange(
    natura_barplot,
    corine_barplot,
    natura_boxplot,
    corine_boxplot,
    ncol = 2
)



par(mfrow = c(3, 4))

for (habitat_type in colnames(fractions_natura)) {
    plot_histogram(env_data_natura[,habitat_type],
                   title = sprintf("%s", habitat_type))
}


par(mfrow = c(5, 4))

for (habitat_type in colnames(fractions_corine)) {
    plot_histogram(env_data_corine[,habitat_type],
                   title = sprintf("%s", habitat_type))
}

par(mfrow = c(3, 2))
plot_histogram(env_data_natura$Temperature,
               title = "Temperature (April and May average, celcius)")
plot_histogram(env_data_natura$Rainfall,
               title = "Rainfall (April and May average, mm)") 
plot_histogram(env_data_natura$PatchDensity,
               title = "Natura patch density within buffer areas")
plot_histogram(env_data_corine$PatchDensity,
               title = "Corine patch density within buffer areas")
plot_histogram(env_data_natura$Effort,
               title = "Transect length (meters)")







# WITHIN DATASET CORRELATIONS
natura_correlation_columns <- c(colnames(fractions_natura),
                                "Temperature",
                                "Rainfall",
                                "PatchDensity")
corine_correlation_columns <- c(colnames(fractions_corine),
                                "Temperature",
                                "Rainfall",
                                "PatchDensity")
pheatmap(cor(env_data_natura[,natura_correlation_columns]), 
         display_numbers = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
pheatmap(cor(env_data_corine[,corine_correlation_columns]), 
         display_numbers = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)



# ABUNDANCE DISTRIBUTIONS
abundance_df <- as.data.frame(abundance)
abundance_long <- stack(abundance_df)
colnames(abundance_long) <- c("Abundance", "Species")

abundance_transects_df <- aggregate(abundance_df, 
                                    list(Transect = spatiotemporal_context$Transect),
                                    FUN = sum)
abundance_transects_df$Transect <- NULL
abundance_transects_long <- stack(abundance_transects_df)
colnames(abundance_transects_long) <- c("Abundance", "Species")

abundance_years_df <- aggregate(abundance_df, 
                                list(Year = spatiotemporal_context$Year),
                                FUN = sum)
abundance_years_df$Year <- NULL
abundance_years_long <- stack(abundance_years_df)
colnames(abundance_years_long) <- c("Abundance", "Species")

total <- ggplot(abundance_long, aes(x = Abundance, color = Species)) +
    geom_freqpoly(binwidth = 5, size = 1) +
    scale_colour_manual(values = rainbow(ncol(abundance))) +
    theme_minimal() +
    labs(title = "Species abundance across samples",
         x = "Abundance", y = "Count") +
    theme(legend.position = "none")
transect <- ggplot(abundance_transects_long, aes(x = Abundance, color = Species)) +
    geom_freqpoly(binwidth = 5, size = 1) +
    scale_colour_manual(values = rainbow(ncol(abundance))) +
    theme_minimal() +
    labs(title = "Species abundance across transects",
         x = "Abundance", y = "Count") +
    theme(legend.position = "none")
year <- ggplot(abundance_years_long, aes(x = Abundance, color = Species)) +
    geom_freqpoly(binwidth = 5, size = 1) +
    scale_colour_manual(values = rainbow(ncol(abundance))) +
    theme_minimal() +
    labs(title = "Species abundance across years",
         x = "Abundance", y = "Count") +
    theme(legend.position = "none")

grid.arrange(
    total, transect, year
)





# PREVALENCE DISTRIBUTIONS

# Prevalence over samples
sample_prevalence <- colSums(occurrence) / nrow(occurrence)

# Prevalence over transects
number_of_transects <- length(unique(spatiotemporal_context$Transect))
transect_prevalence <- sapply(colnames(occurrence), function(species) {
    indices_where_species_is_present <- which(occurrence[, species] == 1)
    transects_where_species_is_present <- unique(spatiotemporal_context$Transect[indices_where_species_is_present])
    transect_prevalence <- length(transects_where_species_is_present) / number_of_transects
    return(transect_prevalence)
})

# Prevalence over years
number_of_years <- length(unique(spatiotemporal_context$Year))
year_prevalence <- sapply(colnames(occurrence), function(species) {
    indices_where_species_is_present <- which(occurrence[, species] == 1)
    years_where_species_is_present <- unique(spatiotemporal_context$Year[indices_where_species_is_present])
    year_prevalence <- length(years_where_species_is_present) / number_of_years
    return(year_prevalence)
})

species_order <- names(sort(sample_prevalence, decreasing = TRUE))
combined_data <- rbind(
    sample = sample_prevalence[species_order],
    transect = transect_prevalence[species_order],
    year = year_prevalence[species_order]
)


old_par <- par(no.readonly = TRUE) 
par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))

par(mfrow = c(1, 3))

barplot(sort(sample_prevalence),
        horiz = TRUE,
        las = 2,
        main = "Sample prevalence")
barplot(sort(transect_prevalence),
        horiz = TRUE,
        las = 2,
        main = "Transect prevalence")
barplot(sort(year_prevalence),
        horiz = TRUE,
        las = 2,
        main = "Year prevalence")

par(old_par)

pheatmap(t(combined_data),
         cluster_rows = FALSE,
         cluster_cols = FALSE)






# SPATIAL AND TEMPORAL PLOTS 


# Average  over transect
natura_transect_averages <- aggregate(env_data_natura, 
                                      by = list(Transect = spatiotemporal_context$Transect),
                                      FUN = mean, na.rm = TRUE)
natura_transect_averages <- merge(natura_transect_averages, 
                                  transect_coordinates, 
                                  by.x = "Transect", 
                                  by.y = "Transect",
                                  all.x = TRUE)
corine_transect_averages <- aggregate(env_data_corine, 
                                      by = list(Transect = spatiotemporal_context$Transect),
                                      FUN = mean, na.rm = TRUE)
corine_transect_averages <- merge(corine_transect_averages, 
                                  transect_coordinates, 
                                  by.x = "Transect", 
                                  by.y = "Transect",
                                  all.x = TRUE)
abundance_transect_averages_by_species <- aggregate(abundance_df,
                                                    by = list(Transect = spatiotemporal_context$Transect),
                                                    FUN = mean)
rownames(abundance_transect_averages_by_species) <- abundance_transect_averages_by_species$Transect
abundance_transect_averages_by_species$Transect <- NULL
abundance_transect_averages <- data.frame(Abundance = rowMeans(abundance_transect_averages_by_species),
                                          Transect = rownames(abundance_transect_averages_by_species))
abundance_transect_averages <- merge(abundance_transect_averages, 
                                     transect_coordinates, 
                                     by.x = "Transect", 
                                     by.y = "Transect",
                                     all.x = TRUE)
transect_species_richness <- c()
for (transect in natura_transect_averages$Transect) {
    row_indices_transect <- which(spatiotemporal_context$Transect == transect)
    occurrences_for_transect <- occurrence[row_indices_transect, , drop=FALSE]
    samples_species_richnesses <- rowSums(occurrences_for_transect)
    species_richness <- mean(samples_species_richnesses)
    transect_species_richness <- c(transect_species_richness, species_richness)
}
transect_species_richness <- data.frame(richness = transect_species_richness,
                                        x = natura_transect_averages$x,
                                        y = natura_transect_averages$y)



# Spatial plots

plot_spatially <- function(data, variable, title, text = FALSE) {
    min_point_size <- 0.5  
    max_point_size <- 3 
    number_of_sizes <- 100
    scaled_values <- scale_to_range(data[, variable], min_point_size, max_point_size)
    plot(data[, c("x")], data[, c("y")], 
         col = "black", 
         pch = 1,
         cex = scaled_values,
         xlab = "X Coordinate", 
         ylab = "Y Coordinate",
         main = title)
    if (text) {
        text(x = data[, c("x")],
             y = data[, c("y")],
             labels = round(data[, variable], 1),  
             pos = 3, 
             col = "black", 
             cex = 0.5)
    }
}

par(mfrow = c(2, 4))
plot_spatially(transect_species_richness, "richness", title = "Average species richness", text=TRUE)
plot_spatially(abundance_transect_averages, "Abundance", title = "Average abundance", text=TRUE)
plot_spatially(natura_transect_averages, "Temperature", title = "Average temperature", text=TRUE)
plot_spatially(natura_transect_averages, "Rainfall", title = "Average rainfall", text=TRUE)
plot_spatially(natura_transect_averages, "PatchDensity", title = "Natura Patch Density", text=TRUE)
plot_spatially(corine_transect_averages, "PatchDensity", title = "Corine Patch Density", text=TRUE)
plot_spatially(natura_transect_averages, "Effort", title = "Transect length", text=TRUE)


# Spatial plots for habitat types
par(mfrow = c(3, 4))
for (type in colnames(fractions_natura)) {
    plot_spatially(natura_transect_averages, type, title = type, text=TRUE)
}

par(mfrow = c(4, 5))
for (type in colnames(fractions_corine)) {
    plot_spatially(corine_transect_averages, type, title = type, text=TRUE)
}



# Average environmental variables over year
natura_year_averages <- aggregate(env_data_natura, 
                                  by = list(Year = spatiotemporal_context$Year),
                                  FUN = mean, na.rm = TRUE)

average_year_abundance_by_species <- aggregate(abundance_df,
                                  by = list(Year = spatiotemporal_context$Year),
                                  FUN = mean)
average_year_abundance <- average_year_abundance_by_species
rownames(average_year_abundance) <- average_year_abundance$Year
average_year_abundance$Year <- NULL
average_year_abundance <- rowMeans(average_year_abundance)
year_species_richness <- c()
transects_sampled <- c()
for (year in natura_year_averages$Year) {
    row_indices_year <- which(spatiotemporal_context$Year == year)
    occurrences_for_year <- occurrence[row_indices_year, , drop=FALSE]
    samples_species_richnesses <- rowSums(occurrences_for_year)
    species_richness <- mean(samples_species_richnesses)
    year_species_richness <- c(year_species_richness, species_richness)
    transects_sampled <- c(transects_sampled, length(row_indices_year))
}


# Temporal plots

par(mfrow=c(4,1))
    
plot(natura_year_averages$Year, 
     natura_year_averages$Temperature, 
     type="l", 
     col="black", 
     lwd = 3,
     ylim=c(-5,5), 
     ylab="Average temperature (°C)",
     xlab="Year")

plot(natura_year_averages$Year, 
     natura_year_averages$Rainfall, 
     type="l", 
     col="black", 
     lwd = 3,
     ylim=c(0,50), 
     ylab="Average rainfall (mm)",
     xlab="Year")

plot(natura_year_averages$Year, 
     year_species_richness, 
     type="l", 
     col="black", 
     lwd = 3,
     ylim=c(0,30), 
     ylab="Average species richness", 
     xlab="Year")

plot(natura_year_averages$Year, 
     transects_sampled, 
     type="l", 
     col="black", 
     lwd = 3,
     ylim=c(0,30), 
     ylab="Transects sampled", 
     xlab="Year")



# SPECIES COMPOSITION
natura_heatmap <- pheatmap(fractions_natura)
row_order_of_natura_heatmap <- natura_heatmap$tree_row$order
transect_order_of_natura_heatmap <- rownames(fractions_natura)[row_order_of_natura_heatmap]
abundance_ordered <- abundance_transect_averages_by_species[transect_order_of_natura_heatmap,]
natura_abundance_heatmap <- pheatmap(log(abundance_ordered + 1),
                                     cluster_rows = FALSE)

corine_heatmap <- pheatmap(fractions_corine)
row_order_of_corine_heatmap <- corine_heatmap$tree_row$order
transect_order_of_corine_heatmap <- rownames(fractions_corine)[row_order_of_corine_heatmap]
abundance_ordered <- abundance_transect_averages_by_species[transect_order_of_corine_heatmap,]
corine_abundance_heatmap <- pheatmap(log(abundance_ordered + 1),
                                     cluster_rows = FALSE)






# SITE FIDELITY
# SITE FIDELITY
species_site_fidelities <- c()
species_prevalences <- c()
for (species in colnames(occurrence)) {
    species_prevalence <- 0
    site_fidelities <- c()
    for (transect in unique(spatiotemporal_context$Transect)) {
        transect_visit_indices <- spatiotemporal_context$Transect == transect
        transect_observations <- occurrence[transect_visit_indices, , drop = FALSE]
        times_visited <- sum(transect_visit_indices)
        times_observed <- sum(transect_observations[, species] == 1)
        if (times_observed == 0) {
            next;
        }
        species_prevalence <- species_prevalence + 1 
        site_fidelity <- times_observed / times_visited
        site_fidelities <- c(site_fidelities, site_fidelity)
    }
    average_site_fidelity <- mean(site_fidelities)
    species_site_fidelities <- c(species_site_fidelities, average_site_fidelity)
    species_prevalences <- c(species_prevalences, species_prevalence)
}
names(species_site_fidelities) <- colnames(occurrence)
names(species_prevalences) <- colnames(occurrence)

sorted_indices <- order(species_site_fidelities)
species_site_fidelities <- species_site_fidelities[sorted_indices]
species_prevalences <- species_prevalences[sorted_indices]

color_palette <- colorRampPalette(c("red", "blue"))(length(species_prevalences))
colors <- color_palette[rank(species_prevalences)]

old_par <- par(no.readonly = TRUE) 
par(mfrow = c(1, 2))
par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))
barplot(species_site_fidelities,
        horiz = TRUE,
        las = 2,
        col = colors,
        main = "Species site fidelity (colored by prevalence)",
        xlab = "Site Fidelity")
plot(species_prevalences, species_site_fidelities)
par(old_par)











# CLOSE PDF
dev.off()
