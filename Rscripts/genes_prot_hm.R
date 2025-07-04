library(tidyverse) 
library(pheatmap)  
nucleotide_csv_path <- "path/to/csvfile.csv" 
protein_csv_path <- "path/to/csvfile.csv"   
output_folder <- "path/to/folder" 

cat("Checking paths and output folder...\n")
cat("Nucleotide path:", nucleotide_csv_path, " - Exists:", file.exists(nucleotide_csv_path), "\n")
cat("Protein path:", protein_csv_path, " - Exists:", file.exists(protein_csv_path), "\n")

if (!dir.exists(output_folder)) {
  cat("The output folder does not exist. Trying to create it:", output_folder, "\n")
  dir_created <- dir.create(output_folder, recursive = TRUE, showWarnings = TRUE)
  if (dir_created) {
    cat("Output folder successfully created.\n")
  } else {
    stop("Unable to create output folder. Check permissions: ", output_folder)
  }
} else {
  cat("Output folder already exists:", output_folder, "\n")
}

cat("Starting processing of CSV files...\n")

# --- 2. Processing for the Nucleotide CSV ---
cat("\nProcessing nucleotide data...\n")

tryCatch({
  # Read the nucleotide CSV file - CHANGED delim = ";" to sep = ";"
  # If you have other packages loaded with a read_csv function having different parameters,
  # you may need to specify readr::read_csv to force using the tidyverse version.
  df_nt_long <- read.csv(nucleotide_csv_path, sep = ";", stringsAsFactors = FALSE) # Using base read.csv, with sep
  
  # Add a check for column names, as read.csv may not infer them well,
  # and also for identity conversion
  cat("Nucleotide file read. Dimensions:", dim(df_nt_long), "\n")
  
  # Check if the necessary columns exist
  if (!all(c("gene", "reference", "identity") %in% names(df_nt_long))) {
    stop("The nucleotide file does not contain columns 'gene', 'reference', and 'identity'. Found column names: ", paste(names(df_nt_long), collapse = ", "))
  }
  
  # Convert the 'identity' column to numeric, handling the comma as a decimal separator
  df_nt_long <- df_nt_long %>%
    mutate(identity = as.numeric(gsub(",", ".", as.character(identity))))
  cat("Nucleotide 'identity' column converted to numeric.\n")
  
  # Pivot from "long" to "wide" format (matrix)
  df_nt_matrix_data <- df_nt_long %>%
    pivot_wider(names_from = reference, values_from = identity) %>%
    column_to_rownames("gene")
  
  cat("Nucleotide data pivoted. Dataframe dimensions:", dim(df_nt_matrix_data), "\n")
  
  # Convert the dataframe to a numeric matrix for pheatmap
  nt_matrix <- as.matrix(df_nt_matrix_data)
  cat("Nucleotide matrix created. Dimensions:", dim(nt_matrix), "\n")
  
  # Check if the matrix is empty or contains only NA
  if (nrow(nt_matrix) == 0 || ncol(nt_matrix) == 0 || all(is.na(nt_matrix))) {
    stop("The nucleotide matrix is empty or contains only NA values. Unable to create heatmap.")
  }
  
  # Create and save the Nucleotide Heatmap
  output_nt_path <- file.path(output_folder, "nucleotide_heatmap.png")
  cat("Attempting to save nucleotide heatmap to:", output_nt_path, "\n")
  png(output_nt_path, width=900, height=800, res=150)
  pheatmap(nt_matrix,
           main = "Nucleotide Identity Heatmap of Key Genes",
           color = colorRampPalette(c("white", "lightblue"))(100),
           cluster_rows = TRUE,    
           cluster_cols = TRUE,    
           display_numbers = TRUE, 
           number_format = "%.1f", 
           fontsize_number = 8,    
           na_col = "grey"         
  )
  dev.off() 
  cat("Nucleotide heatmap successfully saved.\n")
  
}, error = function(e) {
  cat("!!! Critical error in processing nucleotide file: ", e$message, "\n")
  if (!is.null(dev.list())) {
    dev.off() 
    cat("Graphics device closed due to error.\n")
  }
})


# --- 3. Processing for the Protein (Amino Acid) CSV ---
cat("\nProcessing amino acid data...\n")

tryCatch({
  # Read the protein CSV file - CHANGED delim = ";" to sep = ";"
  df_aa_long <- read.csv(protein_csv_path, sep = ";", stringsAsFactors = FALSE) # Using base read.csv, with sep
  
  cat("Amino acid file read. Dimensions:", dim(df_aa_long), "\n")
  
  # Check if the necessary columns exist
  if (!all(c("protein", "reference", "identity") %in% names(df_aa_long))) {
    stop("The protein file does not contain columns 'protein', 'reference', and 'identity'. Found column names: ", paste(names(df_aa_long), collapse = ", "))
  }
  
  # Convert the 'identity' column to numeric, handling the comma as a decimal separator
  df_aa_long <- df_aa_long %>%
    mutate(identity = as.numeric(gsub(",", ".", as.character(identity))))
  
  cat("Amino acid 'identity' column converted to numeric.\n")
  
  # Pivot from "long" to "wide" format (matrix)
  df_aa_matrix_data <- df_aa_long %>%
    pivot_wider(names_from = reference, values_from = identity) %>%
    column_to_rownames("protein")
  
  cat("Amino acid data pivoted. Dataframe dimensions:", dim(df_aa_matrix_data), "\n")
  
  # Convert the dataframe to a numeric matrix for pheatmap
  aa_matrix <- as.matrix(df_aa_matrix_data)
  cat("Amino acid matrix created. Dimensions:", dim(aa_matrix), "\n")
  
  # Check if the matrix is empty or contains only NA
  if (nrow(aa_matrix) == 0 || ncol(aa_matrix) == 0 || all(is.na(aa_matrix))) {
    stop("The amino acid matrix is empty or contains only NA values. Unable to create heatmap.")
  }
  
  # Create and save the Protein Heatmap
  output_aa_path <- file.path(output_folder, "aminoacid_heatmap.png")
  cat("Attempting to save amino acid heatmap to:", output_aa_path, "\n")
  png(output_aa_path, width=900, height=800, res=150)
  pheatmap(aa_matrix,
           main = "Amino Acid Identity Heatmap of Key Genes",
           color = colorRampPalette(c("white", "lightblue"))(100),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.1f",
           fontsize_number = 8,
           na_col = "grey"
  )
  dev.off() 
  cat("Amino acid heatmap successfully saved.\n")
  
}, error = function(e) {
  cat("!!! Critical error in processing protein file: ", e$message, "\n")
  if (!is.null(dev.list())) {
    dev.off() 
    cat("Graphics device closed due to error.\n")
  }
})

cat("\nProcess completed. Check the folder '", output_folder, "' for your heatmaps.\n", sep = "")
