# FUNCTIONS | TRANSECTS AND SPECIES DATA

get_transect_species_data <- function(observations, type) {
    species_list <- unique(observations$Species)
    transect_species_data_list <- list()
    for (transect in unique(observations$Transect)) {
        observations_for_transect <- observations[observations$Transect == transect,]
        data_for_transect <- get_species_data_for_transect(observations_for_transect,
                                                           transect,
                                                           type,
                                                           species_list)
        transect_species_data_list <- append(transect_species_data_list,
                                             list(data_for_transect))
    }
    transect_species_data <- do.call(rbind, transect_species_data_list)
    rownames(transect_species_data) <- transect_species_data$Transect
    transect_species_data$Transect <- NULL
    
    return(transect_species_data)
}


get_species_data_for_transect <- function(data_for_transect, 
                                          transect, 
                                          type, 
                                          species_list) {
    
    if (type == "fraction") {
        data_for_transect <- get_species_fractions_data(data_for_transect, 
                                                        species_list,
                                                        transect)
    } else if (type == "diversity") {
        data_for_transect <- get_species_diversity_data(data_for_transect, 
                                                        transect)
    } else if (type == "presence") {
        data_for_transect <- get_species_presence_data(data_for_transect, 
                                                       species_list,
                                                       transect)
    }
    
    return(data_for_transect)
}



get_species_fractions_data <- function(transect_observations, species_list, transect) {
    species_abundances <- aggregate(Abundance ~ Species, 
                                    data = transect_observations, 
                                    FUN = sum)
    abundance_fractions <- frequencies_to_fractions(species_abundances$Abundance)
    transposed_abundances <- as.data.frame(t(species_abundances[,c("Abundance"), drop = FALSE]))
    colnames(transposed_abundances) <- species_abundances$Species
    transposed_abundances <- rbind(transposed_abundances, Fraction = abundance_fractions)
    for (species in species_list) {
        if (!species %in% colnames(transposed_abundances)) {
            transposed_abundances[,species] <- 0
        }
    }
    transposed_abundances$Transect <- transect
    species_fractions_data <- transposed_abundances["Fraction",, drop = FALSE]
    
    return(species_fractions_data)
}


get_species_diversity_data <- function(data_for_transect, transect) {
    
    species_table <- xtabs(Abundance ~ Species, data_for_transect)
    genus_table <- xtabs(Abundance ~ Genus, data_for_transect)
    species_in_transect <- unique(data_for_transect$Species)
    richness <- c(length(unique(data_for_transect$Species)))
    richness_genus <- c(length(unique(data_for_transect$Species)))
    simpsons_diversity <-  diversity(species_table, "simpson")
    simpsons_diversity_genus <- diversity(genus_table, "simpson")
    shannons_diversity <- diversity(species_table, "shannon")
    shannons_diversity_genus <- diversity(genus_table, "shannon")
    data_for_transect <- data.frame(Transect = c(transect), 
                                    Richness = richness,
                                    RichnessGenus = richness_genus,
                                    SimpsonsDiversity = simpsons_diversity,
                                    SimpsonsDiversityGenus = simpsons_diversity_genus,
                                    ShannonsDiversity = shannons_diversity,
                                    ShannonsDiversityGenus = shannons_diversity_genus)
    
    return(data_for_transect)
}


get_species_presence_data <- function(data_for_transect, species_list, transect) {
    presence_data <- data.frame(Transect = c(transect))
    for (species in species_list) {
        if (species %in% unique(data_for_transect$Species)) {
            presence_data[,species] <- 1
        } else {
            presence_data[,species] <- 0
        }
    }
    data_for_transect <- presence_data
    
    return(data_for_transect)
}


# FUNCTIONS | OTHER




plot_all_habitat_counts <- function(habitat_raster, habitat_names) {
    habitat_counts <- freq(habitat_raster)
    habitat_counts_with_names <- merge(habitat_counts, 
                                       habitat_names, 
                                       by.x = "value", 
                                       by.y = "value",
                                       all.x = TRUE)
    par(mai=c(1,2,1,1))
    barplot(habitat_counts_with_names$count,
            horiz = TRUE,
            names.arg = habitat_counts_with_names$name,
            cex.names = 0.8,
            main = "Habitat type counts",
            xlab = "Count",
            las = 1)
}


cluster_data <- function(data_types, transect_data_list, cluster_amount, category) {
    
    for (data_type in data_types) {
        data_on_transects <- transect_data_list[[data_type]]
        data_on_transects$Cluster <- NULL
        if (data_type == "landscapemetrics" || data_type == "diversity") {
            data_on_transects <- scale_between_zero_and_one(data_on_transects)
        }
        clustering <- hclust(dist(data_on_transects), method = "complete")
        cluster_groups <- cutree(tree = as.dendrogram(clustering), k = cluster_amount)
        annotation <- data.frame(cluster_groups)
        rownames(annotation) <- rownames(data_on_transects)
        pheatmap(data_on_transects,
                 annotation_row = annotation,
                 main = sprintf("%s %s on transects, buffer: %f m", 
                                category, 
                                data_type, 
                                buffer_width))
        data_on_transects$Cluster <- cluster_groups
        transect_data_list[[data_type]] <- data_on_transects
    }
    return(transect_data_list)
}


compare_habitats_and_species <- function(habitat_data_types,
                                         transect_habitat_data_list,
                                         transect_species_data_list) {
    for (habitat_data_type in habitat_data_types) {
        habitat_data_on_transects <- transect_habitat_data_list[[habitat_data_type]]
        habitat_data_on_transects$Cluster <- NULL
        if (habitat_data_type == "landscapemetrics") {
            habitat_data_on_transects <- scale_between_zero_and_one(habitat_data_on_transects)
        }
        habitat_heatmap <- pheatmap(habitat_data_on_transects, 
                                    main = sprintf("Comparison: Habitat %s, buffer: %f m", habitat_data_type, buffer_width))
        row_order_of_heatmap_data <- habitat_heatmap$tree_row$order
        transect_names_in_order_of_heatmap_data <- rownames(habitat_data_on_transects)[row_order_of_heatmap_data]
        for (species_data_type in species_data_types) {
            species_data_on_transects <- transect_species_data_list[[species_data_type]]
            species_data_on_transects$Cluster <- NULL
            if (species_data_type == "diversity") {
                species_data_on_transects <- scale_between_zero_and_one(species_data_on_transects)
            }
            species_data_on_transects_ordered <- species_data_on_transects[transect_names_in_order_of_heatmap_data,]
            pheatmap(species_data_on_transects_ordered, 
                     cluster_rows = FALSE,
                     main = sprintf("Comparison: Species %s, buffer: %f m", species_data_type, buffer_width))
        }
    }
}





# EXPLORE RAW DATA

# Check that landscape is suitable for the metrics
#check_landscape(habitats) # Habitats ok!
#print(list_lsm(level = "landscape"), n = 66) # list all available metrics

# Open pdf for plots exploratory analysis
pdf(file.path(dir_results, "exploratory_analysis.pdf"))

# set filename for exploratory analysis text file
file <- file.path(dir_results, "exploratory_analysis.txt")


# LOAD RAW DATA
load(file = file.path(dir_data, "occurrence_raw.RData")) 
load(file = file.path(dir_data, "fractions_natura.RData"))
load(file = file.path(dir_data, "fractions_corine.RData")) 
load(file = file.path(dir_data, "clusters_natura.RData"))
load(file = file.path(dir_data, "clusters_corine.RData")) 
load(file = file.path(dir_data, "env_data_natura_raw.RData")) 
load(file = file.path(dir_data, "env_data_corine_raw.RData")) 
load(file = file.path(dir_data, "corine_diversities.RData"))
load(file = file.path(dir_data, "natura_diversities.RData"))
load(file = file.path(dir_data, "spatiotemporal_context_raw.RData"))



# BASIC NUMBERS ON RAW DATA


# Occurrence
number_of_species <- length(colnames(occurrence))
append_to_file(sprintf("Number of species: %s\n", number_of_species), file)
# Natura habitats
number_of_habitats_natura <- length(unique(natura_raster)[,1])
append_to_file(sprintf("Number of natura habitats: %s\n", number_of_habitats_natura), file)
number_of_known_habitats_within_buffer_natura <- length(colnames(fractions_natura))
append_to_file(sprintf("Number of known natura habitats within buffer: %s\n", 
                       number_of_known_habitats_within_buffer_natura), file)
# Corine habitats
number_of_habitats_corine <- length(unique(corine_raster)[,1])
append_to_file(sprintf("Number of corine habitats: %s\n", number_of_habitats_corine), file)
number_of_known_habitats_within_buffer_corine <- length(colnames(fractions_corine))
append_to_file(sprintf("Number of known corine habitats within buffer: %s\n", 
                       number_of_known_habitats_within_buffer_corine), file)



# BASIC PLOTS FROM RAW DATA

# On how many years was each transect visited 
transect_visits <- aggregate(SampleEvent ~ Transect, data = observations,
                             FUN = function(x) length(unique(x)))
plot(as.factor(transect_visits$Transect), 
     transect_visits$SampleEvent, 
     type="l",
     main="Number of visits for each transect")
hist(transect_visits$SampleEvent, main = "Histogram for number of visits per transect")

# Plot amount of transects that was visited each year
year_transect_amounts <- aggregate(Transect ~ Year, data = observations,
                                   FUN = function(x) length(unique(x)))
plot(year_transect_amounts$Year, year_transect_amounts$Transect, type="h")
hist(year_transect_amounts$Transect, main = "Histogram for number of transects visited each year")

# Plot transect lengths
transect_names <- unique(observations
                         $Transect)
transect_lengths <- c()
for (transect in transect_names) {
    transect_lengths <- c(transect_lengths, mean(observations
                                                 [observations
                                                     $Transect == transect, "Effort"]))
}
transect_efforts <- data.frame(Transect = transect_names, Effort = transect_lengths)
hist(transect_efforts$Effort, main = "Histogram of transect lengths")






# Produce data
transect_natura_data_list <- list()
transect_corine_data_list <- list()
for (habitat_data_type in habitat_data_types) {
    natura_data_on_transects <- get_transect_habitat_data(buffer_width, 
                                                          natura_raster, 
                                                          transects_shp,
                                                          names_natura,
                                                          habitat_data_type)
    transect_natura_data_list[[habitat_data_type]] <- natura_data_on_transects
    corine_data_on_transects <- get_transect_habitat_data(buffer_width,
                                                          corine_raster,
                                                          transects_shp,
                                                          names_corine,
                                                          habitat_data_type)
    transect_corine_data_list[[habitat_data_type]] <- corine_data_on_transects
}
transect_species_data_list <- list()
for (species_data_type in species_data_types) {
    species_data_on_transects <- get_transect_species_data(observations,
                                                           species_data_type)
    transect_species_data_list[[species_data_type]] <- species_data_on_transects
}


# Cluster data
# TO DO: fix the rest of the code so cluster groups can be saved as column
# Linnustoklusterit lintulajeille (columns)
transect_natura_data_list <- cluster_data(habitat_data_types, 
                                          transect_natura_data_list,
                                          10,
                                          "Habitat")
transect_corine_data_list <- cluster_data(habitat_data_types, 
                                          transect_corine_data_list,
                                          10,
                                          "Habitat")
transect_species_data_list <- cluster_data(species_data_types,
                                           transect_species_data_list,
                                           5,
                                           "Species")
number_of_transects <- length(rownames(transect_natura_data_list[["fraction"]]))
ordered_transects <- rownames(transect_natura_data_list[["fraction"]])[order(transect_natura_data_list[["fraction"]]$Cluster)]
natura_cluster <- transect_natura_data_list[["fraction"]][ordered_transects,]$Cluster
corine_cluster <- transect_corine_data_list[["fraction"]][ordered_transects,]$Cluster
cluster_plot_data <- data.frame(transect = factor(c(ordered_transects, ordered_transects), 
                                                  levels = ordered_transects),
                                cluster = c(natura_cluster, corine_cluster),
                                data = c(rep("natura", times = number_of_transects),
                                         rep("corine", times = number_of_transects)))

ggplot(cluster_plot_data, aes(x = transect, y = cluster, color = data)) +
    geom_point(alpha = 0.7) +
    labs(title = "Transect natura and corine clusters",
         x = "Transect",
         y = "Cluster") +
    theme_minimal() +
    scale_color_manual(values = c("natura" = "blue", "corine" = "red")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))




# Plot transects sampling years and species richness for each transect species combination
all_years <- rev(sort(unique(observations$Year)))
all_transects <- unique(observations$Transect)
sampling_years <- data.frame(years = all_years)
species_richnesses <- data.frame(years = all_years)
transect_lengths <- data.frame(years = all_years)
transect_clusters_natura <- data.frame(years = all_years)
transect_clusters_corine <- data.frame(years = all_years)
for (transect in unique(observations$Transect)) {
    transect_years <- unique(observations[observations$Transect == transect, ]$Year)
    transect_sampling_years <- c()
    transect_species_richness <- c()
    transect_effort <- c()
    transect_cluster_natura <- c()
    transect_cluster_corine <- c()
    for (year in all_years) {
        transect_sampling_years <- c(transect_sampling_years, year %in% transect_years)
        transect_species_richness <- c(transect_species_richness, 
                                       length(unique(observations[observations$Transect == transect & observations$Year == year,]$Species)))
        if (year %in% transect_years) {
            transect_effort <- c(transect_effort, 
                                 transect_efforts[transect_efforts$Transect == transect,]$Effort)
            natura_fractions <- transect_natura_data_list[["fraction"]]
            corine_fractions <- transect_corine_data_list[["fraction"]]
            transect_cluster_natura <- c(transect_cluster_natura,
                                         natura_fractions[transect,]$Cluster)
            transect_cluster_corine <- c(transect_cluster_corine,
                                         corine_fractions[transect,]$Cluster)
        } else {
            transect_effort <- c(transect_effort, 0)
            transect_cluster_natura <- c(transect_cluster_natura, 0)
            transect_cluster_corine <- c(transect_cluster_corine, 0)
        }
        
    }
    sampling_years[,transect] <- transect_sampling_years
    species_richnesses[,transect] <- transect_species_richness
    transect_lengths[,transect] <- transect_effort
    transect_clusters_natura[,transect] <- transect_cluster_natura
    transect_clusters_corine[,transect] <- transect_cluster_corine
}
numeric_sampling_years <- as.data.frame(sapply(sampling_years[-1], as.numeric)) # Exclude the 'years' column
rownames(numeric_sampling_years) <- sampling_years$years
species_richnesses_matrix <- as.matrix(species_richnesses[-1]) # Exclude the first column (years)
rownames(species_richnesses_matrix) <- species_richnesses$years
transect_lengths_matrix <- as.matrix(transect_lengths[-1])
rownames(transect_lengths_matrix) <- transect_lengths$years
transect_clusters_matrix <- as.matrix(transect_clusters_natura[-1])
rownames(transect_clusters_matrix) <- transect_clusters_natura$years
transect_clusters_corine_matrix <- as.matrix(transect_clusters_corine[-1])
rownames(transect_clusters_corine_matrix) <- transect_clusters_corine$years
species_richnesses_matrix <- species_richnesses_matrix[, sort(colnames(numeric_sampling_years))]
transect_lengths_matrix <- transect_lengths_matrix[, sort(colnames(numeric_sampling_years))]
transect_clusters_matrix <- transect_clusters_matrix[, sort(colnames(numeric_sampling_years))]
transect_clusters_corine_matrix <- transect_clusters_corine_matrix[, sort(colnames(numeric_sampling_years))]
pheatmap(numeric_sampling_years,
         cluster_rows = FALSE,
         cluster_columns = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Sampling years for each transect")
pheatmap(species_richnesses_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Species richness for each sample")
pheatmap(transect_lengths_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Transect lengths for each sample")
pheatmap(transect_clusters_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Natura type clusters for each sample")
pheatmap(transect_clusters_corine_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Corine land cover clusters for each sample")


# Plot habitat counts
plot_all_habitat_counts(natura_raster, names_natura)
plot_all_habitat_counts(corine_raster, names_corine)



compare_habitats_and_species(habitat_data_types, 
                             transect_natura_data_list, 
                             transect_species_data_list)
compare_habitats_and_species(habitat_data_types, 
                             transect_corine_data_list, 
                             transect_species_data_list)




# CORRELATIONS FOR CORINE AND NATURA RASTERS

# Calculate correlations between corine and natura data
cor(clusters_natura, clusters_corine)
cor(natura_diversities$SimpsonsDiversity, 
    corine_diversities$SimpsonsDiversity)
cor(natura_diversities$ShannonsDiversity, 
    corine_diversities$ShannonsDiversity)
cor(natura_diversities$PatchDensity, 
    corine_diversities$PatchDensity)
cor(natura_diversities$ScaledRichness, 
    corine_diversities$ScaledRichness)


# Does natura data predict corine and vice versa?
# Do they predict their own clusters? 
lm_data <- data.frame(natura_cluster = clusters_natura,
                      corine_cluster = clusters_corine,
                      fractions_natura,
                      fractions_corine)
formula_natura_frac_to_natura_clus <- as.formula(sprintf("natura_cluster ~ %s", 
                                                         paste(colnames(fractions_natura), collapse = "+")))
formula_corine_frac_to_natura_clus <- as.formula(sprintf("natura_cluster ~ %s", 
                                                         paste(colnames(fractions_corine), collapse = "+")))
formula_corine_frac_to_corine_clus <- as.formula(sprintf("corine_cluster ~ %s", 
                                                         paste(colnames(fractions_corine), collapse = "+")))
formula_natura_frac_to_corine_clus <- as.formula(sprintf("corine_cluster ~ %s", 
                                                         paste(colnames(fractions_natura), collapse = "+")))
model_natura_frac_to_natura_clus <- lm(formula_natura_frac_to_natura_clus, data = lm_data)
model_corine_frac_to_natura_clus <- lm(formula_corine_frac_to_natura_clus, data = lm_data)
model_corine_frac_to_corine_clus <- lm(formula_corine_frac_to_corine_clus, data = lm_data)
model_natura_frac_to_corine_clus <- lm(formula_natura_frac_to_corine_clus, data = lm_data)
model_corine_clus_to_natura_clus <- lm(natura_cluster ~ corine_cluster, data = lm_data)
model_natura_clus_to_corine_clus <- lm(corine_cluster ~ natura_cluster, data = lm_data)
AIC(model_natura_frac_to_natura_clus, model_corine_frac_to_natura_clus)
AIC(model_corine_frac_to_corine_clus, model_natura_frac_to_corine_clus)
AIC(model_corine_clus_to_natura_clus, model_natura_clus_to_corine_clus)


# How well do corine and natura predict occurrences for single species?
lm_data_occurrence <- data.frame(occurrence,
                                 env_data_natura,
                                 env_data_corine,
                                 transect = spatiotemporal_context$Transect,
                                 year = spatiotemporal_context$Year)
models_natura_probit <- list()
models_natura_logit <- list()
models_corine_probit <- list()
models_corine_logit <- list()
for (i in 1:length(colnames(occurrence))) {
    species <- colnames(occurrence)[i]
    print(sprintf("Species %s, %s/%s", species, i, length(colnames(occurrence))))
    formula_natura <- as.formula(sprintf("%s ~ %s", 
                                         species,
                                         paste(colnames(fractions_natura), collapse = "+")))
    formula_corine <- as.formula(sprintf("%s ~ %s", 
                                         species,
                                         paste(colnames(fractions_corine), collapse = "+")))
    models_natura_probit[[i]] <- glm(formula_natura, data = lm_data_occurrence, family = binomial(link = "probit"))
    models_natura_logit[[i]] <- glm(formula_natura, data = lm_data_occurrence, family = binomial(link = "logit"))
    models_corine_probit[[i]] <- glm(formula_corine, data = lm_data_occurrence, family = binomial(link = "probit"))
    models_corine_logit[[i]] <- glm(formula_corine, data = lm_data_occurrence, family = binomial(link = "logit"))
}

wins_corine <- 0
wins_natura <- 0
wins_probit <- 0
wins_logit <- 0
for (i in 1:length(colnames(occurrence))) {
    model_natura_probit <- models_natura_probit[[i]]
    model_natura_logit <- models_natura_logit[[i]]
    model_corine_probit <- models_corine_probit[[i]]
    model_corine_logit <- models_corine_logit[[i]]
    corine_natura_probit <- AIC(model_natura_probit, model_corine_probit)
    corine_natura_logit <- AIC(model_natura_logit, model_corine_logit)
    logit_probit_natura <- AIC(model_natura_logit, model_natura_probit)
    logit_probit_corine <- AIC(model_corine_logit, model_corine_probit)
    if (corine_natura_probit["model_natura_probit",]$AIC < corine_natura_probit["model_corine_probit",]$AIC) {
        wins_natura <- wins_natura + 1 
    } else {
        wins_corine <- wins_corine + 1
    }
    if (corine_natura_logit["model_natura_logit",]$AIC < corine_natura_logit["model_corine_logit",]$AIC) {
        wins_natura <- wins_natura + 1 
    } else {
        wins_corine <- wins_corine + 1
    }
    if (logit_probit_natura["model_natura_probit",]$AIC < logit_probit_natura["model_natura_logit",]$AIC) {
        wins_probit <- wins_probit + 1
    } else {
        wins_logit <- wins_logit + 1
    }
    if (logit_probit_corine["model_corine_probit",]$AIC < logit_probit_corine["model_corine_logit",]$AIC) {
        wins_probit <- wins_probit + 1
    } else {
        wins_logit <- wins_logit + 1
    }
}


for (i in 1:length(models_natura_probit)) {
    print(summary(models_natura_probit[[i]]))
}



# DO SAME FOR SELECTED DATA
load(file = file.path(dir_data, "occurrence.RData")) 
load(file = file.path(dir_data, "env_data_natura.RData")) 
load(file = file.path(dir_data, "env_data_corine.RData")) 



# How well do corine and natura predict occurrences for single species?
lm_data_occurrence <- data.frame(occurrence,
                                 env_data_natura,
                                 env_data_corine)
natura_types_in_selected <- c()
for (type in colnames(fractions_natura)) {
    if (type %in% colnames(env_data_natura)) {
        natura_types_in_selected <- c(natura_types_in_selected, type)
    }
}
corine_types_in_selected <- c()
for (type in colnames(fractions_corine)) {
    if (type %in% colnames(env_data_corine)) {
        corine_types_in_selected <- c(corine_types_in_selected, type)
    }
}


models_natura_probit <- list()
models_natura_logit <- list()
models_corine_probit <- list()
models_corine_logit <- list()
for (i in 1:length(colnames(occurrence))) {
    species <- colnames(occurrence)[i]
    print(sprintf("Species %s, %s/%s", species, i, length(colnames(occurrence))))
    formula_natura <- as.formula(sprintf("%s ~ %s", 
                                         species,
                                         paste(natura_types_in_selected, collapse = "+")))
    formula_corine <- as.formula(sprintf("%s ~ %s", 
                                         species,
                                         paste(corine_types_in_selected, collapse = "+")))
    models_natura_probit[[i]] <- glm(formula_natura, data = lm_data_occurrence, family = binomial(link = "probit"))
    models_natura_logit[[i]] <- glm(formula_natura, data = lm_data_occurrence, family = binomial(link = "logit"))
    models_corine_probit[[i]] <- glm(formula_corine, data = lm_data_occurrence, family = binomial(link = "probit"))
    models_corine_logit[[i]] <- glm(formula_corine, data = lm_data_occurrence, family = binomial(link = "logit"))
}

wins_corine <- 0
wins_natura <- 0
wins_probit <- 0
wins_logit <- 0
for (i in 1:length(colnames(occurrence))) {
    model_natura_probit <- models_natura_probit[[i]]
    model_natura_logit <- models_natura_logit[[i]]
    model_corine_probit <- models_corine_probit[[i]]
    model_corine_logit <- models_corine_logit[[i]]
    corine_natura_probit <- AIC(model_natura_probit, model_corine_probit)
    corine_natura_logit <- AIC(model_natura_logit, model_corine_logit)
    logit_probit_natura <- AIC(model_natura_logit, model_natura_probit)
    logit_probit_corine <- AIC(model_corine_logit, model_corine_probit)
    if (corine_natura_probit["model_natura_probit",]$AIC < corine_natura_probit["model_corine_probit",]$AIC) {
        wins_natura <- wins_natura + 1 
    } else {
        wins_corine <- wins_corine + 1
    }
    if (corine_natura_logit["model_natura_logit",]$AIC < corine_natura_logit["model_corine_logit",]$AIC) {
        wins_natura <- wins_natura + 1 
    } else {
        wins_corine <- wins_corine + 1
    }
    if (logit_probit_natura["model_natura_probit",]$AIC < logit_probit_natura["model_natura_logit",]$AIC) {
        wins_probit <- wins_probit + 1
    } else {
        wins_logit <- wins_logit + 1
    }
    if (logit_probit_corine["model_corine_probit",]$AIC < logit_probit_corine["model_corine_logit",]$AIC) {
        wins_probit <- wins_probit + 1
    } else {
        wins_logit <- wins_logit + 1
    }
}

wins_corine
wins_natura
wins_probit
wins_logit


for (i in 1:length(models_natura_probit)) {
    print(summary(models_natura_probit[[i]]))
}







# CLOSE PDF
dev.off()
