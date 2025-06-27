# THIS FILE CONTAINS FUNCTIONS THAT ARE EITHER USED IN THE MAIN SCRIPT, 
# OR USED BY MULTIPLE SUBSCRIPTS


# GENERAL FUNCTIONS

frequencies_to_fractions <- function(frequencies) {
    total_count <- sum(frequencies)
    fractions <-  (frequencies / total_count)
    
    return(fractions)
} 

append_to_file <- function(text, file, sep="") {
    cat(text, file = file, sep = sep, append = TRUE)
} 

scale_between_zero_and_one <- function(data) {
    scaled_data <- apply(data, 
                         MARGIN = 2, # Apply to columns
                         FUN = function(X) (X - min(X))/diff(range(X)))
    return(as.data.frame(scaled_data))
}

extract_thinning_value <- function(filename) {
    parts <- strsplit(filename, "_")[[1]] 
    thinning_value <- as.numeric(parts[length(parts) - 2])
    print(thinning_value)
    return(thinning_value)
}

scale_to_range <- function(values, new_min, new_max) {
    old_min <- min(values)
    old_max <- max(values)
    scaled_values <- new_min + (values - old_min) * (new_max - new_min) / (old_max - old_min)
    return(scaled_values)
}


plot_histogram <- function(values, title) {
    hist(values, 
         xlab = sprintf("[%s, %s], mean %s, median %s, sd %s",
                        round(min(values), 3), 
                        round(max(values), 3),
                        round(mean(values), 3), 
                        round(median(values), 3),
                        round(sd(values), 3)),
         main = title)
}


# FUNCTIONS | DATA READING

read_climate_data <- function(folder, type) {
    # Set filename prefix based on given climate data type
    prefix <- ''
    if (type == "temperature") {
        prefix <- "tmon"
    }
    if (type == "rainfall") {
        prefix <- "rrmon"
    }
    # List all files starting with the given prefix
    raster_files <- list.files(folder, pattern = sprintf("^%s", prefix), full.names = TRUE)
    # Read rasters in a list
    rasters <- list()
    for (filepath in raster_files) {
        filename <- basename(filepath)
        # Create regex representing the filename. Capture year and month in parenthesis to extract them
        filename_pattern <- sprintf("^%s_(\\d{4})(\\d{2})\\d{2}\\.tif$", prefix)
        year <- sub(filename_pattern, "\\1", filename)
        month <- sub(filename_pattern, "\\2", filename)
        # Create descriptive name for raster and add to list
        raster_name <- sprintf("%s_%s_%s", type, year, month)
        rasters[raster_name] <- list(rast(filepath))
    }
    return(rasters)
}


# FUNCTIONS | TRANSECTS AND HABITATS DATA

get_transect_habitat_data <- function(buffer_width,
                                      habitats,
                                      transects,
                                      habitat_names,
                                      type) {
    all_habitat_values <- unique(habitats)
    all_habitat_values[,1] <- as.numeric(all_habitat_values[,1])
    colnames(all_habitat_values) <- c("value")
    transect_habitat_data_list <- list()
    for (i in 1:length(transects)) {
        transect <- transects[i,]
        transect_number <- transect$Numero[1]
        data_for_transect <- ''
        buffer_polygon_around_transect <- buffer(transect, width = buffer_width)
        if (type == "fraction" || type == "presence") {
            habitats_in_buffer <- na.omit(extract(habitats, buffer_polygon_around_transect))
            habitats_in_buffer <- droplevels(habitats_in_buffer)
            habitats_in_buffer[,2] <- as.integer(as.character(habitats_in_buffer[,2]))
            data_for_transect <- get_habitat_type_fractions_dataframe(transect_number,
                                                                      habitats_in_buffer,
                                                                      all_habitat_values,
                                                                      habitat_names,
                                                                      type)
        } else if (type == "landscapemetrics") {
            habitats_in_buffer <- crop(habitats, 
                                       buffer_polygon_around_transect,
                                       mask = TRUE)
            data_for_transect <- get_habitat_landscapemetrics(transect_number,
                                                              habitats_in_buffer)
        }
        transect_habitat_data_list <- append(transect_habitat_data_list,
                                             list(data_for_transect))
    }
    transect_habitat_data <- do.call(rbind, transect_habitat_data_list)
    rownames(transect_habitat_data) <- transect_habitat_data$Transect
    transect_habitat_data$Transect <- NULL
    
    return(transect_habitat_data)
}


get_habitat_type_fractions_dataframe <- function(transect_number,
                                                 habitat_cell_values_in_buffer,
                                                 all_habitat_values,
                                                 habitat_names,
                                                 type) {
    habitat_type_fractions <- calculate_habitat_type_fractions(transect_number,
                                                               habitat_cell_values_in_buffer,
                                                               all_habitat_values)
    habitat_type_fractions$Transect <- NULL     # Remove transect column to prevent column type change to string during transpose
    if (type == "presence") {
        habitat_type_fractions$Value <- habitat_type_fractions$Value > 0
    }
    habitat_type_fractions_transposed <- t(habitat_type_fractions)
    colnames(habitat_type_fractions_transposed) <- habitat_type_fractions$Habitat
    habitat_type_fractions_transposed <- as.data.frame(habitat_type_fractions_transposed["Value",, drop = FALSE])
    indices_of_colname_values_in_habitat_names <- match(colnames(habitat_type_fractions_transposed), habitat_names$value)
    colnames(habitat_type_fractions_transposed) <- habitat_names$name[indices_of_colname_values_in_habitat_names]
    habitat_type_fractions_transposed$Transect <- transect_number
    
    return(habitat_type_fractions_transposed)
}



calculate_habitat_type_fractions <- function(transect_number, 
                                             habitat_cell_values_in_buffer, 
                                             habitat_types) {
    habitat_type_frequencies <- table(habitat_cell_values_in_buffer)
    habitat_type_frequencies <- as.data.frame(habitat_type_frequencies, stringsAsFactors = FALSE)
    colnames(habitat_type_frequencies) <- c("Transect", "Habitat", "Value")
    habitat_type_frequencies$Transect <- transect_number
    habitat_type_frequencies <- add_value_zero_for_missing_habitat_types(habitat_type_frequencies, 
                                                                         habitat_types, 
                                                                         transect_number)
    habitat_type_frequencies$Habitat <- as.double(habitat_type_frequencies$Habitat)
    habitat_type_frequencies$Value <- as.double(habitat_type_frequencies$Value)
    habitat_type_frequencies$Value <-  frequencies_to_fractions(habitat_type_frequencies$Value)
    
    return(habitat_type_frequencies)
}



add_value_zero_for_missing_habitat_types <- function(habitat_data, habitat_types, transect_number) {
    for (habitat_type in habitat_types$value) {
        if (!as.numeric(habitat_type) %in% habitat_data$Habitat) {
            habitat_data[nrow(habitat_data) + 1,] <- c(transect_number, habitat_type, 0)
        }
    }
    
    return(habitat_data)
    
} 


get_habitat_landscapemetrics <- function(transect_number, habitats) {
    #patch_richness <- lsm_l_pr(habitats)
    patch_richness_per_area <- lsm_l_prd(habitats)
    simpsons_diversity <- lsm_l_sidi(habitats)
    shannons_diversity <- lsm_l_shdi(habitats)
    shannons_evenness <- lsm_l_shei(habitats)
    patch_density <- lsm_l_pd(habitats)
    largest_patch_index <- lsm_l_lpi(habitats)
    edge_density <- lsm_l_ed(habitats)
    transect_landscapemetrics <- data.frame(Transect = transect_number,
                                            #Richness = patch_richness$value,
                                            ScaledRichness = patch_richness_per_area$value,
                                            SimpsonsDiversity = simpsons_diversity$value,
                                            ShannonsDiversity = shannons_diversity$value,
                                            ShannonsEvenness = shannons_evenness$value,
                                            PatchDensity = patch_density$value,
                                            LargetPatchIndex = largest_patch_index$value,
                                            EdgeDensity = edge_density$value)
    
    return(transect_landscapemetrics)
}


