library(tidyverse) # Per read_csv (da readr), pivot_wider, mutate, filter, column_to_rownames
library(pheatmap)  # Per generare le heatmap

# --- 1. Imposta i percorsi dei tuoi file CSV e della cartella di output ---
# *** IMPORTANTE: MODIFICA QUESTI PERCORSI CON QUELLI CORRETTI SUL TUO COMPUTER ***
nucleotidico_csv_path <- "C:/Users/silvi/Desktop/grafici_tesi/heatmap/nuc.csv" 
proteico_csv_path <- "C:/Users/silvi/Desktop/grafici_tesi/heatmap/prot.csv"   
output_folder <- "C:/Users/silvi/Desktop/grafici_tesi/heatmap" # Sembra che tu voglia salvarle nella stessa cartella

# --- DEBUG: Verifica l'esistenza dei percorsi e la creazione della cartella di output ---
cat("Verifica percorsi e output folder...\n")
cat("Percorso nucleotidico:", nucleotidico_csv_path, " - Esiste:", file.exists(nucleotidico_csv_path), "\n")
cat("Percorso proteico:", proteico_csv_path, " - Esiste:", file.exists(proteico_csv_path), "\n")

if (!dir.exists(output_folder)) {
  cat("La cartella di output non esiste. Provando a crearla:", output_folder, "\n")
  dir_created <- dir.create(output_folder, recursive = TRUE, showWarnings = TRUE)
  if (dir_created) {
    cat("Cartella di output creata con successo.\n")
  } else {
    stop("Impossibile creare la cartella di output. Controlla i permessi: ", output_folder)
  }
} else {
  cat("La cartella di output esiste già:", output_folder, "\n")
}

cat("Inizio elaborazione dei file CSV...\n")

# --- 2. Processo per il CSV Nucleotidico ---
cat("\nElaborazione dati nucleotidici...\n")

tryCatch({
  # Leggi il file CSV nucleotidico - CAMBIATO delim = ";" con sep = ";"
  # Se hai caricato altri pacchetti che hanno una funzione read_csv con lo stesso nome
  # ma con parametri diversi, potresti dover specificare readr::read_csv
  # per forzare l'uso della versione di tidyverse.
  df_nt_long <- read.csv(nucleotidico_csv_path, sep = ";", stringsAsFactors = FALSE) # Usiamo read.csv base, con sep
  
  # Aggiungiamo anche il controllo dei nomi delle colonne, perché read.csv potrebbe non inferire bene
  # e per la conversione dell'identità
  # E anche col_types = cols() non è un parametro di read.csv, quindi lo tolgo
  
  cat("File nucleotidico letto. Dimensioni:", dim(df_nt_long), "\n")
  
  # Controlla se le colonne necessarie esistono
  if (!all(c("gene", "reference", "identity") %in% names(df_nt_long))) {
    stop("Il file nucleotidico non contiene le colonne 'gene', 'reference' e 'identity'. Nomi colonne trovati: ", paste(names(df_nt_long), collapse = ", "))
  }
  
  # Trasforma la colonna 'identity' in numerico, gestendo la virgola come separatore decimale
  df_nt_long <- df_nt_long %>%
    mutate(identity = as.numeric(gsub(",", ".", as.character(identity))))
  cat("Colonna 'identity' nucleotidica convertita a numerico.\n")
  
  # Trasforma da formato "lungo" a formato "largo" (matrice)
  df_nt_matrix_data <- df_nt_long %>%
    pivot_wider(names_from = reference, values_from = identity) %>%
    column_to_rownames("gene")
  
  cat("Dati nucleotidici pivottati. Dimensioni del dataframe:", dim(df_nt_matrix_data), "\n")
  
  # Converti il dataframe in una matrice numerica per pheatmap
  nt_matrix <- as.matrix(df_nt_matrix_data)
  cat("Matrice nucleotidica creata. Dimensioni:", dim(nt_matrix), "\n")
  
  # Controlla se la matrice è vuota o contiene solo NA
  if (nrow(nt_matrix) == 0 || ncol(nt_matrix) == 0 || all(is.na(nt_matrix))) {
    stop("La matrice nucleotidica è vuota o contiene solo valori NA. Impossibile creare la heatmap.")
  }
  
  # Crea e salva la Heatmap Nucleotidica
  output_nt_path <- file.path(output_folder, "heatmap_nucleotidica.png")
  cat("Tentativo di salvare heatmap nucleotidica in:", output_nt_path, "\n")
  png(output_nt_path, width=900, height=800, res=150)
  pheatmap(nt_matrix,
           main = "Heatmap identità nucleotidica dei geni chiave",
           color = colorRampPalette(c("white", "lightblue"))(100),
           cluster_rows = TRUE,    
           cluster_cols = TRUE,    
           display_numbers = TRUE, 
           number_format = "%.1f", 
           fontsize_number = 8,    
           na_col = "grey"         
  )
  dev.off() 
  cat("Heatmap nucleotidica salvata con successo.\n")
  
}, error = function(e) {
  cat("!!! Errore critico nell'elaborazione del file nucleotidico: ", e$message, "\n")
  if (!is.null(dev.list())) {
    dev.off() 
    cat("Dispositivo grafico chiuso a causa dell'errore.\n")
  }
})


# --- 3. Processo per il CSV Proteico (Aminoacidico) ---
cat("\nElaborazione dati aminoacidici...\n")

tryCatch({
  # Leggi il file CSV proteico - CAMBIATO delim = ";" con sep = ";"
  df_aa_long <- read.csv(proteico_csv_path, sep = ";", stringsAsFactors = FALSE) # Usiamo read.csv base, con sep
  
  cat("File aminoacidico letto. Dimensioni:", dim(df_aa_long), "\n")
  
  # Controlla se le colonne necessarie esistono
  if (!all(c("protein", "reference", "identity") %in% names(df_aa_long))) {
    stop("Il file proteico non contiene le colonne 'protein', 'reference' e 'identity'. Nomi colonne trovati: ", paste(names(df_aa_long), collapse = ", "))
  }
  
  # Trasforma la colonna 'identity' in numerico, gestendo la virgola come separatore decimale
  df_aa_long <- df_aa_long %>%
    mutate(identity = as.numeric(gsub(",", ".", as.character(identity))))
  
  cat("Colonna 'identity' aminoacidica convertita a numerico.\n")
  
  # Trasforma da formato "lungo" a formato "largo" (matrice)
  df_aa_matrix_data <- df_aa_long %>%
    pivot_wider(names_from = reference, values_from = identity) %>%
    column_to_rownames("protein")
  
  cat("Dati aminoacidici pivottati. Dimensioni del dataframe:", dim(df_aa_matrix_data), "\n")
  
  # Converti il dataframe in una matrice numerica per pheatmap
  aa_matrix <- as.matrix(df_aa_matrix_data)
  cat("Matrice aminoacidica creata. Dimensioni:", dim(aa_matrix), "\n")
  
  # Controlla se la matrice è vuota o contiene solo NA
  if (nrow(aa_matrix) == 0 || ncol(aa_matrix) == 0 || all(is.na(aa_matrix))) {
    stop("La matrice aminoacidica è vuota o contiene solo valori NA. Impossibile creare la heatmap.")
  }
  
  # Crea e salva la Heatmap Proteica
  output_aa_path <- file.path(output_folder, "heatmap_aminoacidica.png")
  cat("Tentativo di salvare heatmap aminoacidica in:", output_aa_path, "\n")
  png(output_aa_path, width=900, height=800, res=150)
  pheatmap(aa_matrix,
           main = "Heatmap identità aminoacidica dei geni chiave",
           color = colorRampPalette(c("white", "lightblue"))(100),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.1f",
           fontsize_number = 8,
           na_col = "grey"
  )
  dev.off() 
  cat("Heatmap aminoacidica salvata con successo.\n")
  
}, error = function(e) {
  cat("!!! Errore critico nell'elaborazione del file proteico: ", e$message, "\n")
  if (!is.null(dev.list())) {
    dev.off() 
    cat("Dispositivo grafico chiuso a causa dell'errore.\n")
  }
})

cat("\nProcesso completato. Controlla la cartella '", output_folder, "' per le tue heatmap.\n", sep = "")
