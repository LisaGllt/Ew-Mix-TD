source("functions/fun.R")
f_load_libraries_colors()

# 1. Design ----

# 5 ratios and 7 levels per ratio + 1 control = 36 conditions

ratios <- c("E", "F", "G", "H", "I") 
levels <- 1:7 
controle <- "N" 
replicats <- 4 

df_exp <- expand.grid(Ratio = ratios, Niveau = levels, Replicat = 1:replicats) |> 
  mutate(Condition = paste0(Ratio, Niveau))

df_ctrl <- data.frame(
  Ratio = controle,
  Niveau = 0,
  Replicat = rep(1:replicats, times = 1),
  Condition = paste0(controle,"0")
)

df_design_exp <- bind_rows(df_exp, df_ctrl) |> 
  mutate(ID_cond_rep = paste(Condition, Replicat, sep="."))

df_design_lotD <- df_design_exp |> 
  filter(Replicat %in% c(1,2))

df_design_lotE <- df_design_exp |> 
  filter(Replicat %in% c(3,4))

# We need 36 culture boxes for each lot

# 2. Culture boxes ----

numbers <- 1:80
colors <- c("Pink", "Orange")

df_culture <- expand.grid(Number = numbers, Color = colors) |>
  mutate(ID_vdt_culture = paste(Number, Color))

# 3. Randomization ----

set.seed(1212)

Nb_boxes_random <- sample(numbers)
Boxes_lotD <- Nb_boxes_random[1:36]
Boxes_lotE <- Nb_boxes_random[37:72]

df_culture_lotD <- df_culture |> 
  filter(Number %in% Boxes_lotD)
df_culture_lotE <- df_culture |> 
  filter(Number %in% Boxes_lotE)

# For lot D

ID_random <- sample(df_culture_lotD$ID_vdt_culture)

ID_unique_D <- unique(Boxes_lotD)
ID_unique_random <- sample(Boxes_lotD)

df_box2 <- data.frame(
  ID_random.1 = ID_unique_D,
  ID_random.2 = ID_unique_random
)

df_vdt_Dtmp <- data.frame(
  Name_rep.1 = df_design_lotD$ID_cond_rep,
  Name_rep.2 = df_design_lotD$ID_cond_rep,
  ID_random.1 = ID_random,
  Earthworm.1 = rep(1, 72),
  Earthworm.2 = rep(2, 72)
) |> 
  mutate(
    Color.1 = sub(".*(Orange|Pink)$", "\\1", ID_random.1),
    Color.2 = case_when(
      Color.1 == "Orange" ~ "Green",
      Color.1 == "Pink" ~ "Yellow"
      ),
    ID_random.1 = as.numeric(substr(ID_random, 1, 2))
  ) |> 
  left_join(df_box2)

df_vdt_1 <- data.frame(
  Name_rep = df_vdt_Dtmp$Name_rep.1,
  Earthworm = df_vdt_Dtmp$Earthworm.1,
  ID_boxe = df_vdt_Dtmp$ID_random.1,
  Color = df_vdt_Dtmp$Color.1
)

df_vdt_2 <- data.frame(
  Name_rep = df_vdt_Dtmp$Name_rep.2,
  Earthworm = df_vdt_Dtmp$Earthworm.2,
  ID_boxe = df_vdt_Dtmp$ID_random.2,
  Color = df_vdt_Dtmp$Color.2
)

# Fot lot E
ID_random <- sample(df_culture_lotE$ID_vdt_culture)

ID_unique_E <- unique(Boxes_lotE)
ID_unique_random <- sample(Boxes_lotE)

df_box2 <- data.frame(
  ID_random.3 = ID_unique_E,
  ID_random.4 = ID_unique_random
)

df_vdt_Etmp <- data.frame(
  Name_rep.3 = df_design_lotE$ID_cond_rep,
  Name_rep.4 = df_design_lotE$ID_cond_rep,
  ID_random.3 = ID_random,
  Earthworm.3 = rep(1, 72),
  Earthworm.4 = rep(2, 72)
) |> 
  mutate(
    Color.3 = sub(".*(Orange|Pink)$", "\\1", ID_random.3),
    Color.4 = case_when(
      Color.3 == "Orange" ~ "Green",
      Color.3 == "Pink" ~ "Yellow"
    ),
    ID_random.3 = as.numeric(substr(ID_random, 1, 2))
  ) |> 
  left_join(df_box2)

df_vdt_3 <- data.frame(
  Name_rep = df_vdt_Etmp$Name_rep.3,
  Earthworm = df_vdt_Etmp$Earthworm.3,
  ID_boxe = df_vdt_Etmp$ID_random.3,
  Color = df_vdt_Etmp$Color.3
)

df_vdt_4 <- data.frame(
  Name_rep = df_vdt_Etmp$Name_rep.4,
  Earthworm = df_vdt_Etmp$Earthworm.4,
  ID_boxe = df_vdt_Etmp$ID_random.4,
  Color = df_vdt_Etmp$Color.4
)

df <- rbind(df_vdt_1, df_vdt_2, df_vdt_3, df_vdt_4)

writexl::write_xlsx(df, path = here::here("protocols/3_Mix_tests/MIX2_Attribution.xlsx"))



