# THIS SCRIPTS FORMATS THE DATA IN SUITABLE FORM FOR HMSC ANALYSES

# FUNCTIONS 

get_transect_weather_data <- function(buffer_width,
                                      weather_rasters,
                                      transects) {
    monthly_weather_for_transects <- data.frame()
    transects <- project(transects, crs(weather_rasters[[1]]))
    for (i in 1:length(transects$Numero)) {
        transect <- transects[i,]
        transect_number <- transects$Numero[i]
        buffer_polygon_around_transect <- buffer(transect, width = buffer_width)
        transect_monthly_mean_weather <- list()
        for (raster_name in names(weather_rasters)) {
            raster <- weather_rasters[[raster_name]]
            raster_name_parts <- strsplit(raster_name, "_")[[1]]
            year <- as.numeric(raster_name_parts[2])
            month <- as.numeric(raster_name_parts[3])
            weather_values_in_buffer <- na.omit(extract(raster, buffer_polygon_around_transect))
            colnames(weather_values_in_buffer) <- c("id", "value")
            month_mean_weather_for_transect <- mean(weather_values_in_buffer$value)
            transect_monthly_mean_weather[["year"]] <- c(transect_monthly_mean_weather[["year"]], year)
            transect_monthly_mean_weather[["month"]] <-  c(transect_monthly_mean_weather[["month"]], month)
            transect_monthly_mean_weather[["value"]] <-  c(transect_monthly_mean_weather[["value"]],
                                                             month_mean_weather_for_transect)
            
        }
        weather_data_for_transect <- data.frame(transect = transect_number, 
                                                    transect_monthly_mean_weather)
        monthly_weather_for_transects <- rbind(monthly_weather_for_transects, weather_data_for_transect)
        
    }
    return(monthly_weather_for_transects)
}

format_observations_for_hmsc <- function(observations) {
    species_list <- unique(observations$Species)
    observation_data_list <- list()
    for (year in unique(observations$Year)) {
        observations_for_year <- observations[observations$Year == year,]
        transect_species_data_list <- list()
        for (transect in unique(observations_for_year$Transect)) {
            observations_for_transect <- observations_for_year[observations_for_year$Transect == transect,]
            data_for_transect <- get_species_abundance_for_transect(observations_for_transect,
                                                                    transect,
                                                                    species_list)
            transect_species_data_list <- append(transect_species_data_list,
                                                 list(data_for_transect))
        }
        transect_species_data_for_year <- do.call(rbind, transect_species_data_list)
        transect_species_data_for_year$Year <- year
        observation_data_list <- append(observation_data_list,
                                        list(transect_species_data_for_year))
    } 
    observation_data <- do.call(rbind, observation_data_list)
    return(observation_data)
}

get_species_abundance_for_transect <- function(data_for_transect, 
                                               transect, 
                                               species_list) {
    abundance_data <- data.frame(Transect = c(transect))
    for (species in species_list) {
        if (species %in% unique(data_for_transect$Species)) {
            abundance_data[,species] <- sum(data_for_transect[data_for_transect$Species == species, "Abundance"])
        } else {
            abundance_data[,species] <- 0
        }
    }
    return(abundance_data)
}

match_to_birdtree <- function(species_names, column_to_match) {
    # Get the index of each species name in the alternative names df
    # Index is NA if name is not found in the specified column (no alternative names)
    match_indices <- match(species_names, species_alternative_names[[column_to_match]])
    # Get alternative names based on the match
    birdtree_species_names <- species_alternative_names$birdtree[match_indices]
    # Replace NA values with the original species names
    birdtree_species_names[is.na(birdtree_species_names)] <- species_names[is.na(birdtree_species_names)]
    return (birdtree_species_names)
}


get_fractions <- function(buffer_width, raster, transects_shp, names) {
    fractions <- get_transect_habitat_data(buffer_width, 
                                           raster, 
                                           transects_shp,
                                           names,
                                           "fraction")
    fractions <- fractions[,!is.na(colnames(fractions))]
    fractions <- fractions[,sapply(fractions, function(x) var(x, na.rm = TRUE) != 0)]
}

get_transect_data_percentages <- function(transects_shp, raster, data_type_names, buffer_width) {
    percentages <- c()
    for (i in 1:length(transects_shp)) {
        transect <- transects_shp[i,]
        transect_number <- transect$Numero[1]
        buffer_polygon_around_transect <- buffer(transect, width = buffer_width)
        habitats_in_buffer <- extract(raster, buffer_polygon_around_transect)
        data_exists <- !is.na(habitats_in_buffer[,2]) & as.integer(habitats_in_buffer[,2]) %in% data_type_names$value 
        percentage <- sum(data_exists) / length(data_exists)
        percentages <- c(percentages, percentage)
    }
    return(percentages)
}

get_transect_clusters <- function(fractions, cluster_amount) {
    
    clustering <- hclust(dist(fractions), method = "complete")
    cluster_groups <- cutree(tree = as.dendrogram(clustering), k = cluster_amount)

    return(cluster_groups)
}





# SCRIPT STARTS



# READ IN DATA
natura_classification <- read_excel(file.path(dir_data, "Ylalappi_luokitus.xls"))
observations <- read.csv(file.path(dir_data, "observations_preprocessed.csv"), sep = ";")
transects_shp <- vect(file.path(dir_data, "transects_preprocessed.shp"))
species_traits <- read_excel(file.path(dir_data, "BirdTraits21112018.xlsx"))
species_alternative_names <- read.csv(file.path(dir_data, "speciesAlternativeNames.txt"), sep = ";")
taxonomy <- read.tree(file.path(dir_data, "tree.txt")) #TO DO: CHOOSE RANDOMLY FROM LIST OF TREES?
corine_classification <- read_excel(file.path(dir_data, "CorineMaanpeite2018Luokat.xls"))


# GENERAL REFORMATTING OF RAW DATA
names_natura <- data.frame(value = natura_classification$Value, name = natura_classification$NaturaTyyppi)
names_corine <- data.frame(value = corine_classification$Value, name = corine_classification$Level4Suo)
observations$Transect <- as.character(observations$Transect)









# FORMAT COMMUNITY DATA Y 
# Rows are samples. Columns are abundances or occurrences

# First generate abundances for each species.
# Separate abundances for each combination of year and transect
abundance_samples <- format_observations_for_hmsc(observations)

# Rename species with names used in birdtree phylogenies
colnames(abundance_samples) <- match_to_birdtree(colnames(abundance_samples), "observations")

# Get the names of included species for new phylogeny if needed
writeLines(colnames(abundance_samples), file.path(dir_data, "species.txt"))

# Match the writing format of species names to that of birdtree
colnames(abundance_samples) <- gsub(" ", "_", colnames(abundance_samples))

# Give each sample a unique id
sample_id <- as.factor(1:nrow(abundance_samples))
rownames(abundance_samples) <- sample_id

# Create community matrix Y for just abundances
abundance <- abundance_samples
abundance$Transect <- NULL
abundance$Year <- NULL
abundance <- as.matrix(abundance)

# Creato community matrix Y for occurrences
occurrence <- 1 * (abundance > 0)







# FORMAT SPATIOTEMPORAL CONTEXT S

# Calculate spatial coordinates for transects
# First define center points for each transect which acts as their location
transect_center_points <- centroids(transects_shp)
# Get coordinates of center points.
transect_coordinates <- crds(transect_center_points, df = TRUE)
rownames(transect_coordinates) <- transect_center_points$Numero
# Calculate spatiotemporal id
spatiotemporal_id <- sprintf("%s%s", abundance_samples$Year, abundance_samples$Transect)

# Save the spatiotemporal context of each sample  
spatiotemporal_context <- data.frame(Sample = sample_id,
                                     Transect = abundance_samples$Transect,
                                     Year = abundance_samples$Year,
                                     YearTransect = spatiotemporal_id,
                                     x = transect_coordinates[abundance_samples$Transect, "x"],
                                     y = transect_coordinates[abundance_samples$Transect, "y"])










# FORMAT ENVIRONMENTAL DATA X
# Rows: samples, columns: env data of each sample

# Get fractions of each habitat type on each transect
fractions_natura <- get_fractions(buffer_width, natura_raster, transects_shp, names_natura)
fractions_corine <- get_fractions(buffer_width, corine_raster, transects_shp, names_corine)
colnames(fractions_natura) <- make.names(colnames(fractions_natura))
colnames(fractions_corine) <- make.names(colnames(fractions_corine))

# Get transect names
transect_names <- rownames(fractions_natura)

# Calculate PCA from raw data
pca_results_natura <- prcomp(fractions_natura, center = TRUE, scale. = TRUE)
pca_results_corine <- prcomp(fractions_corine, center = TRUE, scale. = TRUE)


# Calculate habitat type clusters for some easy comparisons
clusters_natura <- get_transect_clusters(fractions_natura, 5)
clusters_corine <- get_transect_clusters(fractions_corine, 5)



# Calculate transect lengths 
transect_lengths <- c()
for (transect in transect_names) {
    transect_lengths <- c(transect_lengths, 
                          mean(observations[observations$Transect == transect, "Effort"]))
}
names(transect_lengths) <- transect_names


# Calculate transect temperature values
transect_temperatures <- get_transect_weather_data(buffer_width, temperature_data, transects_shp)
transect_temperatures_spring <- transect_temperatures[transect_temperatures$month %in% c(4,5), ]
transect_spring_mean_temps <- aggregate(value ~ transect + year, 
                                          data = transect_temperatures_spring, 
                                          FUN = mean, 
                                          na.rm = TRUE)
sample_spring_mean_temperatures <- c()
for (sample_number in 1:length(spatiotemporal_context$Sample)) {
    sample_transect <- spatiotemporal_context[sample_number, ]$Transect
    sample_year <- spatiotemporal_context[sample_number, ]$Year
    sample_spring_mean_temp <- transect_spring_mean_temps[transect_spring_mean_temps$transect == sample_transect 
                                                                & transect_spring_mean_temps$year == sample_year, ]$value
    if (length(sample_spring_mean_temp) == 0)  {
        sample_spring_mean_temp <- NA
    }
    sample_spring_mean_temperatures <- c(sample_spring_mean_temperatures, sample_spring_mean_temp)
}

transect_rainfall <- get_transect_weather_data(buffer_width, rainfall_data, transects_shp)
transect_rainfall_spring <- transect_rainfall[transect_rainfall$month %in% c(4, 5), ]
transect_spring_mean_rainfall <- aggregate(value ~ transect + year, 
                                           data = transect_rainfall_spring, 
                                           FUN = mean,
                                           na.rm = TRUE)
sample_spring_mean_rainfalls <- c()
for (sample_number in 1:length(spatiotemporal_context$Sample)) {
    sample_transect <- spatiotemporal_context[sample_number, ]$Transect
    sample_year <- spatiotemporal_context[sample_number, ]$Year
    sample_spring_mean_rainfall <- transect_spring_mean_rainfall[transect_spring_mean_rainfall$transect == sample_transect 
                                                          & transect_spring_mean_rainfall$year == sample_year, ]$value
    if (length(sample_spring_mean_rainfall) == 0)  {
        sample_spring_mean_rainfall <- NA
    }
    sample_spring_mean_rainfalls <- c(sample_spring_mean_rainfalls, sample_spring_mean_rainfall)
}



# Create diversity data
natura_diversities <- get_transect_habitat_data(buffer_width,
                                                natura_raster,
                                                transects_shp,
                                                names_natura,
                                                "landscapemetrics")
corine_diversities <- get_transect_habitat_data(buffer_width,
                                                corine_raster,
                                                transects_shp,
                                                names_corine,
                                                "landscapemetrics")

# Calculate the percentage of transect over which there is existing habitat type data
transect_natura_data_percentages <- get_transect_data_percentages(transects_shp, 
                                                                  natura_raster, 
                                                                  names_natura,
                                                                  buffer_width)
names(transect_natura_data_percentages) <- transect_names
transect_corine_data_percentages <- get_transect_data_percentages(transects_shp, 
                                                                  corine_raster, 
                                                                  names_corine,
                                                                  buffer_width)
names(transect_corine_data_percentages) <- transect_names



# Create environmental data X for natura
env_data_natura <- data.frame(fractions_natura,
                              Effort = transect_lengths,
                              PatchDensity = natura_diversities$PatchDensity,
                              SimpsonsDiversity = natura_diversities$SimpsonsDiversity,
                              ShannonsDiversity = natura_diversities$ShannonsDiversity,
                              ScaledRichness = natura_diversities$ScaledRichness,
                              NaturaPercentage = transect_natura_data_percentages,
                              CorinePercentage = transect_corine_data_percentages,
                              Cluster = clusters_natura)
env_data_natura <- env_data_natura[spatiotemporal_context$Transect, ]
env_data_natura$Temperature <- sample_spring_mean_temperatures
env_data_natura$Rainfall <- sample_spring_mean_rainfalls
env_data_corine <- data.frame(fractions_corine,
                              Effort = transect_lengths,
                              PatchDensity = corine_diversities$PatchDensity,
                              SimpsonsDiversity = corine_diversities$SimpsonsDiversity,
                              ShannonsDiversity = corine_diversities$ShannonsDiversity,
                              ScaledRichness = corine_diversities$ScaledRichness,
                              NaturaPercentage = transect_natura_data_percentages,
                              CorinePercentage = transect_corine_data_percentages,
                              Cluster = clusters_corine)
env_data_corine <- env_data_corine[spatiotemporal_context$Transect, ]
env_data_corine$Temperature <- sample_spring_mean_temperatures
env_data_corine$Rainfall <- sample_spring_mean_rainfalls













# FORMAT TRAIT DATA T
# Rows: species, columns: traits

trait_data <- species_traits

# Write species names to trait data in the format of the taxonomy tree
trait_data$Species <- gsub(" ", "_", trait_data$Sp_Lat)

# Check that all species have traits
observed_species <- colnames(abundance)
for (species in observed_species) {
    found <- species %in% unique(trait_data$Species)
    if (found == FALSE) {
        print(species)
    }
}

# Include only observed species in trait data
trait_data <- as.data.frame(trait_data[trait_data$Species %in% observed_species,])

# Transform mass into logmass
trait_data$LogMass <- log(trait_data$Mass)

# Name rows based on species
rownames(trait_data) <- trait_data$Species


# Add species phylogenetic order as trait

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
                 "Käkilinnut",
                 "Kyyhkylinnut")
groups <- cutree(taxonomy, 12) # Jaa lahkotasolle (manuaalisesti selvitetty)

species_orders <- c()
for (species in rownames(trait_data)) {
    species_order <- groups[species]
    species_order_name <- order_names[species_order]
    species_orders <- c(species_orders, species_order_name)
}
trait_data$Order <- species_orders


# Order species to same order as they are in observation data
trait_data <- trait_data[observed_species,]







# FORMAT PHYLOGENETIC DATA C
phylogeny_data <- taxonomy






# SAVE DATA

# Save raw data
#save(natura_raster, file = file.path(dir_data, "natura_raster.RData"))
#save(corine_raster, file = file.path(dir_data, "corine_raster.RData"))
save(names_natura, file = file.path(dir_data, "names_natura.RData"))
save(names_corine, file = file.path(dir_data, "names_corine.RData"))
save(observations, file = file.path(dir_data, "observations.RData"))
writeVector(transects_shp, 
            file.path(dir_data, "transects_formatted.shp"),
            overwrite = TRUE)
save(species_traits, file = file.path(dir_data, "species_traits.RData"))
save(species_alternative_names, file = file.path(dir_data, "species_alternative_names.RData"))
save(taxonomy, file = file.path(dir_data, "taxonomy.RData"))
save(transect_lengths, file = file.path(dir_data, "transect_lengths.RData"))
#save(temperature_data, file = file.path(dir_data, "temperature_data.RData"))
#save(rainfall_data, file = file.path(dir_data, "rainfall_data.RData"))

# Save formatted data
save(fractions_natura, file = file.path(dir_data, "fractions_natura_raw.RData"))
save(fractions_corine, file = file.path(dir_data, "fractions_corine_raw.RData"))
save(pca_results_natura, file = file.path(dir_data, "pca_results_natura_raw.RData"))
save(pca_results_corine, file = file.path(dir_data, "pca_results_corine_raw.RData"))
save(clusters_natura, file = file.path(dir_data, "clusters_natura.RData"))
save(clusters_corine, file = file.path(dir_data, "clusters_corine.RData"))
save(transect_lengths, file = file.path(dir_data, "transect_lengths.RData"))
save(transect_temperatures, file = file.path(dir_data, "transect_temperatures.RData"))
save(natura_diversities, file = file.path(dir_data, "natura_diversities.RData"))
save(corine_diversities, file = file.path(dir_data, "corine_diversities.RData"))
save(transect_natura_data_percentages, file = file.path(dir_data, "transect_natura_data_percentages.RData"))
save(transect_corine_data_percentages, file = file.path(dir_data, "transect_corine_data_percentages.RData"))
save(abundance_samples, file = file.path(dir_data, "abundance_samples.RData"))
save(transect_coordinates, file = file.path(dir_data, "transect_coordinates.RData"))

# Save HMSC data
save(env_data_natura, file = file.path(dir_data, "env_data_natura_raw.RData"))
save(env_data_corine, file = file.path(dir_data, "env_data_corine_raw.RData"))
save(abundance, file = file.path(dir_data, "abundance_raw.RData"))
save(occurrence, file = file.path(dir_data, "occurrence_raw.RData"))
save(spatiotemporal_context, file = file.path(dir_data, "spatiotemporal_context_raw.RData"))
save(trait_data, file = file.path(dir_data, "trait_data_raw.RData"))
save(phylogeny_data, file = file.path(dir_data, "phylogeny_data_raw.RData"))

