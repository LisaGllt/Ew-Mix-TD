library(here)
source(file = here::here("functions/fun.R"))
f_load_libraries_colors()

df_rescape <- read.csv(here::here("data/Data_Rescape_long.csv")) |> 
  subset(pesticide %in% c("EPX", "IMD"))

df_rescape_stats <-  df_rescape |> 
  group_by(pesticide, trt)|> 
  summarise(
    median = median(dose),
    mean = mean(dose),
    C_max = max(dose),
    C_min = min(dose)
  ) 

df_rescape_stats_long <- df_rescape_stats |> 
  pivot_longer(cols = c(median, mean, C_max, C_min)) |> 
  mutate(
    habitat = "general",
    agri_type = case_when(
      trt == "nt" ~ "OF",
      trt == "tt" ~ "conventionnal"
    ),
    dose = value,
    type_conc = name,
    source = "Rescape project"
  ) |> 
  dplyr::select(pesticide, habitat, agri_type, type_conc, dose, source)

df_litt_ecotox <- read_excel(here::here("data/Data_litt_EPX_IMD.xlsx"), sheet=1)
df_litt_soils <- read_excel(here::here("data/Data_litt_EPX_IMD.xlsx"), sheet=2)

df_litt_soils <- rbind(df_litt_soils, df_rescape_stats_long)

p <- ggplot(
  data = df_rescape,
  aes(
    x = trt,
    y = dose,
    color = pesticide,
    fill = pesticide
  )
)+
  geom_boxplot(alpha=0.2)+
  geom_dotplot(
    binaxis='y', 
    stackdir='center', 
    dotsize=0.7, 
    alpha=0.7
  )+
  facet_wrap(~pesticide, ncol=2, scale="free")+
  scale_color_manual(values = col_molec)+
  scale_fill_manual(values = col_molec)+
  theme_minimal()
p

# | fig-cap: Comparison of soil contamination indicators and ecotoxicological parameters for EPX and IMD. Ecotoxicological parameters are from the PPDB and the soil contamination indicators are from @pelosi_2021 and @froger_2023.
# | label: fig-soilconta
# | width: 3
# | height: 6

df_litt_soils_plot <- df_litt_soils |> 
  dplyr::filter(
    type_conc %in% c("median", "C_max", "Q25", "Q75", "C_min"),
    agri_type  %in% c("conventionnal", "general"))

df_litt_soils_plot$type_conc <- fct_relevel(df_litt_soils_plot$type_conc,"C_max", "Q75", "median", "Q25", "C_min")

sizetitle <- 12
shape_froger <-  22
shape_pelosi <- 23
shape_rescape <- 24
dot_size <-  3

p <- ggplot() + 
  
  # Rescape median
  geom_point(
    data = df_litt_soils_plot,
    aes(
      x=" ", 
      y = dose, 
      colour = type_conc,
      shape=source, 
      fill = type_conc
      ),
    size = dot_size,
    #stroke=1.2
    ) +
  
  facet_wrap(~ pesticide+source) +
  
  scale_color_manual(
    name = "",
    values = c(Nord_aurora, Nord_frost, Nord_polar)
    ) +
  scale_fill_manual(
    name = "",
    values = c(Nord_aurora, Nord_frost, Nord_polar)
    ) +
  
  # scale_shape_manual(
  #   values = c(shape_rescape, shape_pelosi, shape_froger),
  #   #label = c("Froger et al. 2023", "Pelosi et al. 2021", "Rescape project")
  #   ) +
  labs(y = "Dose (ng/g)") +
  scale_y_log10()+
  theme_minimal()

p


