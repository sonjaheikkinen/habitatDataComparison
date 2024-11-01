# THIS SCRIPTS FORMATS THE DATA IN SUITABLE FORM FOR HMSC ANALYSES

# FUNCTIONS | DATA READING

read_climate_data <- function(folder, type) {
    orefix <- ''
    if (type == "temperature") {
        prefix <- "tmon"
    }
    if (type == "rainfall") {
        prefix <- "rrmon"
    }
    raster_files <- list.files(folder, pattern = sprintf("^%s", prefix), full.names = TRUE)
    rasters <- c()
    for (filename in raster_files) {
        rasters <- c(rasters, rast(filename))
    }
    return(rasters)
}

# FUNCTIONS | DATA FORMATTING

get_transect_temperature_data <- function(buffer_width,
                                          temperature_rasters,
                                          transects) {
    transect_temperature_data <- c()
    for (i in 1:length(transects$Numero)) {
        transect <- transects[i,]
        transect_number <- transect$Numero[1]
        buffer_polygon_around_transect <- buffer(transect, width = buffer_width)
        all_temp_values <- c()
        for (raster in temperature_rasters) {
            temps_in_buffer <- na.omit(extract(raster, buffer_polygon_around_transect))
            colnames(temps_in_buffer) <- c("id", "value")
            all_temp_values <- c(all_temp_values, unique(temps_in_buffer$value))
        }
        mean_temperature <- mean(all_temp_values)
        transect_temperature_data <- c(transect_temperature_data, mean_temperature)
        
    }
    return(transect_temperature_data)
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


# SCRIPT STARTS
###################################################################################################################

# READ IN DATA
natura_raster <- rast(file.path(dir_data, "natura_2393.tif"))
natura_classification <- read_excel(file.path(dir_data, "Ylalappi_luokitus.xls"))
observations <- read.csv(file.path(dir_data, "observations_preprocessed.csv"), sep = ";")
transects_shp <- vect(file.path(dir_data, "transects_preprocessed.shp"))
species_traits <- read_excel(file.path(dir_data, "BirdTraits21112018.xlsx"))
species_alternative_names <- read.csv(file.path(dir_data, "speciesAlternativeNames.txt"), sep = ";")
taxonomy <- read.tree(file.path(dir_data, "tree.txt")) #TO DO: CHOOSE RANDOMLY FROM LIST OF TREES?
temperature_data <- read_climate_data(file.path(dir_data, "ilmastoaineisto"), "temperature")
rainfall_data <- read_climate_data(file.path(dir_data, "ilmastoaineisto"), "rainfall")


# GENERAL REFORMATTING OF RAW DATA
names_natura <- data.frame(value = natura_classification$Value, name = natura_classification$NaturaTyyppi)
observations$Transect <- as.character(observations$Transect)



# FORMAT ENVIRONMENTAL DATA X
# Rows: samples, columns: env data of each sample

# Get fractions of each habitat type on each transect
fractions_natura <- get_transect_habitat_data(buffer_width, 
                                              natura_raster, 
                                              transects_shp,
                                              names_natura,
                                              "fraction")
fractions_natura <- fractions_natura[,!is.na(colnames(fractions_natura))]
fractions_natura <- fractions_natura[,sapply(fractions_natura, function(x) var(x, na.rm = TRUE) != 0)]

# Get transect names
transect_names <- rownames(fractions_natura)

# Calculate transect lengths 
transect_lengths <- c()
for (transect in transect_names) {
    transect_lengths <- c(transect_lengths, 
                          mean(observations[observations$Transect == transect, "Effort"]))
}
names(transect_lengths) <- transect_names


# Calculate transect temperature values
transect_temperatures <- get_transect_temperature_data(buffer_width,
                                                       temperature_data,
                                                       transects_shp)
names(transect_temperatures) <- transect_names


# Create diversity data
natura_diversities <- get_transect_habitat_data(buffer_width,
                                                natura_raster,
                                                transects_shp,
                                                names_natura,
                                                "landscapemetrics")


# Create environmental data X for natura
env_data_natura <- data.frame(fractions_natura,
                              Effort = transect_lengths,
                              Temperature = transect_temperatures,
                              Diversity = natura_diversities$PatchDensity)
colnames(env_data_natura) <- make.names(env_data_natura)





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



# FORMAT TRAIT DATA T

# FORMAT PHYLOGENETIC DATA C

# FORMAT SPATIOTEMPORAL CONTEXT S

# ORDER ALL DATA IN SAMPLE ORDER


# SAVE DATA

# Save raw data
save(natura_raster, file = file.path(dir_data, "natura_raster.RData"))
save(names_natura, file = file.path(dir_data, "names_natura.RData"))
save(observations, file = file.path(dir_data, "observations.RData"))
save(transects_shp, file = file.path(dir_data, "transects_shp.RData"))
save(species_traits, file = file.path(dir_data, "species_traits.RData"))
save(species_alternative_names, file = file.path(dir_data, "species_alternative_names.RData"))
save(taxonomy, file = file.path(dir_data, "taxonomy.RData"))
save(temperature_data, file = file.path(dir_data, "temperature_data.RData"))
save(rainfall_data, file = file.path(dir_data, "rainfall_data.RData"))

# Save formatted data
save(fractions_natura, file = file.path(dir_data, "fractions_natura.RData"))
save(transect_lengths, file = file.path(dir_data, "transect_lengths.RData"))
save(transcet_temperatures, file = file.path(dir_data, "transect_temperatures.RData"))
save(natura_diversities, file = file.path(dir_data, "natura_diversities.RData"))
save(abundance_samples, file = file.path(dir_data, "abundance_samples.RData"))

# Save HMSC data
save(env_data_natura, file = file.path(dir_data, "env_data_natura.RData"))
