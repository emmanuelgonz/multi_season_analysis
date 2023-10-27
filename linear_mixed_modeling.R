# Load the necessary packages
library(asreml)
library(readr)
library(dplyr)
library(tibble)  # for rownames_to_column function
library(rptR)
library(heritability)

calculate_heritabilities <- function(file_paths) {
  # Initialize a list to store the heritabilities and repeatabilities for all files
  all_heritabilities <- list()
  
  # Loop over each file path
  for (file_path in file_paths) {
    tryCatch({
      # Open the CSV file
      data <- read_csv(file_path)
      
      # Convert the variables to factors
      data$treatment <- as.factor(data$treatment)
      data$genotype <- as.factor(data$genotype)
      # data$rep <- as.factor(data$rep)
      # data$range <- as.factor(data$range)
      # data$column <- as.factor(data$column)
      
      # List of phenotypes
      phenotypes <- colnames(data)[!(colnames(data) %in% c("date", "days_after_planting", "plot",
                                                           "treatment", "genotype", "rep", "range", 
                                                           "column", "index", "year", "min_x", "min_y",
                                                           "min_z", "max_x", "max_y", "max_z", "Unnamed: 0", "season",
                                                           "experiment", "field", "plant_name"))]
      print(phenotypes)
      
      # Remove rows where any phenotype is NA
      data <- data[complete.cases(data[phenotypes]), ]
      
      # Convert all phenotypes to numeric
      data[phenotypes] <- lapply(data[phenotypes], function(x) as.numeric(as.character(x)))
      
      # Initialize a list to store the models
      models <- list()
      
      # Loop over each phenotype
      for (i in phenotypes) {
        # Specify the model
        model <- asreml(as.formula(paste(i, "~ treatment")),
                        random = ~ genotype + genotype:treatment + rep/treatment + range/rep:treatment + column/rep:treatment,
                        data = data)
        # ai.sing = TRUE)
        
        # Get the BLUPs for genotype
        genotype_blups <- predict(model)$random$genotype
        
        # Create a file name with the phenotype name and specify the directory
        file_name <- paste("blups/genotype_blups_", i, ".csv", sep = "")
        
        # Save the BLUPs to a file
        write.csv(genotype_blups, file = file_name)
        
        # Store the model in the list
        models[[i]] <- model
      }
      
      
      # # Loop over each phenotype
      # 
      # for (i in phenotypes) {
      #   model <- tryCatch({
      #     # Specify the model
      #     asreml(as.formula(paste(i, "~ treatment")),
      #            random = ~ genotype + genotype:treatment + rep/treatment + range/rep:treatment + column/rep:treatment,
      #            data = data)
      #   }, error = function(e) {
      #     # If the first model fails, fit a different model without rep/treatment in range and column
      #     asreml(as.formula(paste(i, "~ treatment")),
      #            random = ~ genotype + genotype:treatment + rep/treatment + range/treatment + column/treatment,
      #            data = data,
      #            ai.sing = TRUE)
      #   })
      #   
      #   # Store the model in the list
      #   models[[i]] <- model
      # }
      
      # Initialize a list to store the heritability and repeatability estimates
      heritabilities <- list()
      
      # Calculate the number of unique replicates (R) and treatments (L)
      R <- length(unique(data$rep))
      L <- length(unique(data$treatment))
      
      # Loop over each model
      for (i in names(models)) {
        # Extract the variance components
        varcomp <- summary(models[[i]])$varcomp
        
        # Calculate the genetic variance (Vg)
        Vg <- varcomp["genotype", "component"]
        
        # Calculate the genotype by treatment interaction variance (Vge)
        Vge <- varcomp["genotype:treatment", "component"]
        
        # Calculate the residual variance (Ve)
        Ve <- varcomp["units", "component"]
        
        # Calculate the total phenotypic variance (Vp)
        Vp <- Vg + Vge/L + Ve/(R*L)
        
        # Calculate broad-sense heritability (H^2)
        H2 <- Vg / Vp
        
        # Calculate repeatability (r) using the formula from the image
        r <- Vg / (Vg + Vge/L + Ve/R)
        
        # Store the heritability and repeatability estimates in the list
        heritabilities[[i]] <- c(H2, r)
      }
      
      # Transform heritabilities into a data frame and convert row names to a column
      heritabilities_df <- as.data.frame(do.call(rbind, heritabilities)) %>%
        rownames_to_column(var = "trait")
      
      # Extract date from file path and keep only YYYY-MM-DD format
      date_string <- gsub(".+/(\\d{4}-\\d{2}-\\d{2}).+", "\\1", file_path)
      
      # Add date column to heritabilities_df
      heritabilities_df$date <- date_string
      
      # Rename columns in heritabilities_df
      colnames(heritabilities_df) <- c("trait", "broad_sense_heritability", "repeatability", "date")
      
      # Append heritabilities_df to all_heritabilities
      all_heritabilities[[file_path]] <- heritabilities_df
      
    }, error = function(e) {
      message(paste("Error processing file:", file_path, "\n", e))
    })
  }
  
  # Combine all data frames in the list into one data frame
  all_heritabilities <- do.call(rbind, all_heritabilities)
  
  return(all_heritabilities)
}

# Check if 'heritability' directory exists, if not create it
if (!dir.exists("heritability")) {
  dir.create("heritability")
}

# Check if 'blups' directory exists, if not create it
if (!dir.exists("blups")) {
  dir.create("blups")
}

# Specify your directory path here
dir_path <- "./phenotype_data_by_date"

# Search for all CSV files in the given path
# file_paths <- list.files(path = dir_path, pattern = "\\_individual_plant.csv$", full.names = TRUE)
file_paths <- list.files(path = dir_path, pattern = "\\_plot_average.csv$", full.names = TRUE)

# Call the function with your list of CSV file paths
all_heritabilities <- calculate_heritabilities(file_paths)

# Create output file path with 'heritability/' folder
output_file_path <- "heritability/all_heritabilities.csv"

# Write all_heritabilities to CSV file
write.csv(all_heritabilities, file = output_file_path, row.names = FALSE)  # add row.names = FALSE to exclude row names from output file