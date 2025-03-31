# LOAD DATA
load(file = file.path(dir_data, "abundance_raw.RData"))
load(file = file.path(dir_data, "occurrence_raw.RData"))
load(file = file.path(dir_data, "spatiotemporal_context_raw.RData"))
load(file = file.path(dir_data, "env_data_natura_raw.RData"))
load(file = file.path(dir_data, "fractions_natura_raw.RData"))
load(file = file.path(dir_data, "fractions_corine_raw.RData"))
load(file = file.path(dir_data, "env_data_corine_raw.RData"))
load(file = file.path(dir_data, "trait_data_raw.RData"))
load(file = file.path(dir_data, "phylogeny_data_raw.RData"))


# SELECT TRANSECTS
# Remove transect if not a big enough portion of the transect is covered by natura data
transect_selection <- env_data_natura$NaturaPercentage >= 0.9
abundance <- abundance[transect_selection, ]
occurrence <- occurrence[transect_selection, ]
spatiotemporal_context <- spatiotemporal_context[transect_selection, ]
env_data_natura <- env_data_natura[transect_selection, ]
env_data_corine <- env_data_corine[transect_selection, ]
fractions_natura <- fractions_natura[unique(spatiotemporal_context$Transect),]
fractions_corine <- fractions_corine[unique(spatiotemporal_context$Transect),]



# SELECT SPECIES
# Select only more prevalent species
species_list <- c()
samples <- c()
transects <- c()
prevalences <- c()
for (species in colnames(occurrence)) {
    species_occurrence <- occurrence[,species]
    number_of_transects_where_species_has_occurrence <- sum(tapply(species_occurrence, 
                                                                   spatiotemporal_context$Transect, 
                                                                   function(x) any(x == 1)))
    number_of_samples_species_is_found_in <- sum(species_occurrence)
    prevalence <- number_of_samples_species_is_found_in / nrow(occurrence)
    if (number_of_transects_where_species_has_occurrence < 4) {
        species_list <- c(species_list, species)
        samples <- c(samples, number_of_samples_species_is_found_in)
        transects <- c(transects, number_of_transects_where_species_has_occurrence)
        prevalences <- c(prevalences, prevalence)
        occurrence <- occurrence[, colnames(occurrence) != species]
        abundance <- abundance[, colnames(abundance) != species]
    }
}
values <- data.frame(species = species_list,
                     samples = samples,
                     transects = transects,
                     prevalences = prevalences)
print(values)


# SELECT ENVIRONMENTAL DATA
natura_types <- colnames(fractions_natura)
corine_types <- colnames(fractions_corine)
for (type in natura_types) {
    #zero_percentage <- sum(env_data_natura[,type] == 0) / length(env_data_natura[,type])
    #print(sprintf("%s zero percentage: %s", type, zero_percentage))
    #if (zero_percentage >= 0.95) {
    #    env_data_natura[,type] <- NULL
    #    fractions_natura[,type] <- NULL
    #}
    non_zeros <- sum(fractions_natura[,type] != 0)
    print(sprintf("%s non-zeros: %s", type, non_zeros))
    if (non_zeros < 3) {
            env_data_natura[,type] <- NULL
            fractions_natura[,type] <- NULL
    }
}
for (type in corine_types) {
    #zero_percentage <- sum(env_data_corine[,type] == 0) / length(env_data_corine[,type])
    #print(sprintf("%s zero percentage: %s", type, zero_percentage))
    #if (zero_percentage >= 0.95) {
    #    env_data_corine[,type] <- NULL
    #    fractions_corine[,type] <- NULL
    #}
    non_zeros <- sum(fractions_corine[,type] != 0)
    print(sprintf("%s non-zeros: %s", type, non_zeros))
    if (non_zeros < 3) {
        env_data_corine[,type] <- NULL
        fractions_corine[,type] <- NULL
    }
}

# CALCULATE NEW PCA FROM SELECTED DATA
pca_results_natura <- prcomp(fractions_natura, center = TRUE, scale. = TRUE)
pca_results_corine <- prcomp(fractions_corine, center = TRUE, scale. = TRUE)


# SELECT TRAIT DATA
trait_data <- trait_data[colnames(occurrence), ]
trait_data <- trait_data[,c("Feeding"), drop = FALSE]

# SELECT PHYLOGENY DATA
species_to_keep <- colnames(occurrence)
species_to_remove <- setdiff(phylogeny_data$tip.label, species_to_keep)
phylogeny_data <- drop.tip(phylogeny_data, species_to_remove)


# SAVE DATA
save(env_data_natura, file = file.path(dir_data, "env_data_natura.RData"))
save(env_data_corine, file = file.path(dir_data, "env_data_corine.RData"))
save(abundance, file = file.path(dir_data, "abundance.RData"))
save(occurrence, file = file.path(dir_data, "occurrence.RData"))
save(spatiotemporal_context, file = file.path(dir_data, "spatiotemporal_context.RData"))
save(trait_data, file = file.path(dir_data, "trait_data.RData"))
save(phylogeny_data, file = file.path(dir_data, "phylogeny_data.RData"))
save(fractions_natura, file = file.path(dir_data, "fractions_natura.RData"))
save(fractions_corine, file = file.path(dir_data, "fractions_corine.RData"))
save(pca_results_natura, file = file.path(dir_data, "pca_results_natura.RData"))
save(pca_results_corine, file = file.path(dir_data, "pca_results_corine.RData"))


