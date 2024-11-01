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






# EXPLORE RAW DATA

# Check that landscape is suitable for the metrics
#check_landscape(habitats) # Habitats ok!
#print(list_lsm(level = "landscape"), n = 66) # list all available metrics

# Open pdf for plots exploratory analysis
pdf(file.path(dir_results, "exploratory_analysis.pdf"))

# set filename for exploratory analysis text file
file <- file.path(dir_results, "exploratory_analysis.txt")


# LOAD RAW DATA
load(file = file.path(dir_data, "occurrence_raw.RData")) # occurrence
load(file = file.path(dir_data, "natura_raster.RData")) # natura_raster
load(file = file.path(dir_data, "fractions_natura.RData")) # fractions_natura



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
for (habitat_data_type in habitat_data_types) {
    natura_data_on_transects <- get_transect_habitat_data(buffer_width, 
                                                          natura_raster, 
                                                          transects_shp,
                                                          names_natura,
                                                          habitat_data_type)
    transect_natura_data_list[[habitat_data_type]] <- natura_data_on_transects
}
transect_species_data_list <- list()
for (species_data_type in species_data_types) {
    species_data_on_transects <- get_transect_species_data(observations,
                                                           species_data_type)
    transect_species_data_list[[species_data_type]] <- species_data_on_transects
}




# Close pdf
dev.off()
