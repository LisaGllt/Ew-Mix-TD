library(tidyverse)
library(stringr)
library(writexl)

# Files path
folder_path <- "data/MIX_coc_size"

# List of all .csv files
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Empty data.frame initialisation
all_eggs <- data.frame()

# Initialisation of IDs
egg_id <- 1

# For all files
for (file in csv_files) {
  
  # File name without the extension
  filename <- tools::file_path_sans_ext(basename(file))
  
  #Ratio (A), Line (X) et Nb_rep (Y)
  Ratio <- str_sub(filename, 1, 1)
  Line <- as.numeric(str_sub(filename, 2, 2))
  Nb_rep <- as.numeric(str_sub(filename, 3, 3))
  
  # read file
  df <- read.csv(file)
  
  # Length correspond alternatively to Height then Width
  # Separation in 2 columns
  if (nrow(df) %% 2 != 0) {
    warning(paste("Odd measures in", file, "- Ignored files"))
    next
  }
  
  coc_height <- df$Length[seq(1, nrow(df), by = 2)]
  coc_width  <- df$Length[seq(2, nrow(df), by = 2)]
  
  # Temporary data frame
  temp <- data.frame(
    ID = egg_id:(egg_id + length(coc_height) - 1),
    Ratio = Ratio,
    Line = Line,
    Nb_rep = Nb_rep,
    Coc_height = coc_height,
    Coc_width = coc_width
  )
  
  # ID update
  egg_id <- egg_id + length(coc_height)
  
  # Addition to principal data frame
  all_eggs <- bind_rows(all_eggs, temp)
}

# Voir le résultat
print(all_eggs)

write_csv(all_eggs, "data/Data_cocoon_sizes.csv")
write_xlsx(all_eggs, path = "data/Data_cocoon_sizes.xlsx")

