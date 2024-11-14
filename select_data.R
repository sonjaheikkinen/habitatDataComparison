# LOAD DATA
load(file = file.path(dir_data, "abundance_raw.RData"))
load(file = file.path(dir_data, "occurrence_raw.RData"))
load(file = file.path(dir_data, "spatiotemporal_context_raw.RData"))
load(file = file.path(dir_data, "env_data_natura_raw.RData"))
load(file = file.path(dir_data, "env_data_corine_raw.RData"))
load(file = file.path(dir_data, "trait_data_raw.RData"))


# SELECT YEARS
year_selection <- spatiotemporal_context$Year >= 2006
abundance <- abundance[year_selection, ]
occurrence <- occurrence[year_selection, ]
spatiotemporal_context <- spatiotemporal_context[year_selection, ]
env_data_natura <- env_data_natura[year_selection, ]
env_data_corine <- env_data_corine[year_selection, ]

# SELECT TRANSECTS
# Remove transect if not a big enough portion of the transect is covered by natura data
transect_selection <- env_data_natura$NaturaPercentage >= 0.9
abundance <- abundance[transect_selection, ]
occurrence <- occurrence[transect_selection, ]
spatiotemporal_context <- spatiotemporal_context[transect_selection, ]
env_data_natura <- env_data_natura[transect_selection, ]
env_data_corine <- env_data_corine[transect_selection, ]


# SELECT SPECIES
# Select only more prevalent species
for (species in colnames(occurrence)) {
    species_occurrence <- occurrence[,species]
    species_prevalence <- sum(species_occurrence == 1) / length(species_occurrence)
    print(sprintf("Prevalence %s: %s", species, species_prevalence))
    if (species_prevalence < 0.05) {
        occurrence <- occurrence[, colnames(occurrence) != species]
    }
}


# SELECT ENVIRONMENTAL DATA
not_habitat_types <- c("Effort", "Diversity", "Temperature", "NaturaPercentage", "CorinePercentage", "Cluster")
natura_types <- setdiff(colnames(env_data_natura), not_habitat_types)
corine_types <- setdiff(colnames(env_data_corine), not_habitat_types)
for (type in natura_types) {
    zero_percentage <- sum(env_data_natura[,type] == 0) / length(env_data_natura[,type])
    print(sprintf("%s zero percentage: %s", type, zero_percentage))
    if (zero_percentage > 0.95) {
        env_data_natura[,type] <- NULL
    }
}
for (type in corine_types) {
    zero_percentage <- sum(env_data_corine[,type] == 0) / length(env_data_corine[,type])
    print(sprintf("%s zero percentage: %s", type, zero_percentage))
    if (zero_percentage > 0.95) {
        env_data_corine[,type] <- NULL
    }
}
for (type in not_habitat_types) {
    env_data_natura[,type] <- NULL
    env_data_corine[,type] <- NULL
}

# SELECT TRAIT DATA
trait_data <- trait_data[, c("Feeding", "LogMass")]
trait_data <- trait_data[colnames(occurrence), ]


# SAVE DATA
save(env_data_natura, file = file.path(dir_data, "env_data_natura.RData"))
save(env_data_corine, file = file.path(dir_data, "env_data_corine.RData"))
save(abundance, file = file.path(dir_data, "abundance.RData"))
save(occurrence, file = file.path(dir_data, "occurrence.RData"))
save(spatiotemporal_context, file = file.path(dir_data, "spatiotemporal_context.RData"))
save(trait_data, file = file.path(dir_data, "trait_data.RData"))
save(phylogeny_data, file = file.path(dir_data, "phylogeny_data.RData"))


