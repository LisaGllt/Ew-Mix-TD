library(here)
source(file = here::here("functions/fun.R"))
library(directlabels)
library(geomtextpath)
library(ggrepel)
library(ggtext)

f_load_libraries_colors()

path_fig <- here::here("fig/")


# df_Area <- data.frame(
#   Area = c("Europe", "France", "Germany", "Austria", "Belgium", "Bulgaria", "Croatia", 
#            "Denmark", "Spain", "Estonia", "Finland", "Greece", "Hungary", "Ireland", 
#            "Italia", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands (Kingdom of the)",
#            "Poland", "Portugal", "Czechia", "Romania", "Slovakia", "Slovenia", "Sueden"),
#   Area_french = c("Europe", "France", "Allemagne", "Autriche", "Belgique", "Bulgarie", "Croatie", 
#                   "Danemark", "Espagne", "Estonie", "Finlande", "Grèce", "Hongrie", "Irelande", 
#                   "Italia", "Lettonie", "Lituanie", "Luxembourg", "Malte", "Pays-Bas",
#                   "Pologne", "Portugal", "République tchèque", "Roumanie", "Slovaquie", "Slovénie", "Suède")
# )

df_Area <- data.frame(
  Area = c("Europe", "France", "Germany", "Austria", "Belgium", 
           "Denmark", "Spain", "Luxembourg", "Finland", "Sueden",
           "Italia", "Netherlands (Kingdom of the)", "Ireland",
          "Portugal", "Greece", "United Kingdom of Great Britain and Northern Ireland"),
  Area_french = c("Europe", "France", "Allemagne", "Autriche", "Belgique", 
                  "Danemark", "Espagne", "Luxembourg", "Finlande", "Suède",
                  "Italia", "Pays-Bas", "Irelande",
                  "Portugal", "Grèce", "Royaume-Uni")
)

df_data_FAO <- read_excel(
  here::here("data/Data_pesticide_use.xlsx"), 
  sheet = 2
) |> 
  filter(Area %in% df_Area$Area) |> 
  left_join(df_Area, by = "Area") |> 
  mutate(
    Year = as.numeric(Year),
    Value = as.numeric(Value),
    Color_plot = case_when(
      Area == "France" ~ Nord_frost[2],
      Area == "Europe" ~ Nord_frost[4],
      .default = darken(Nord_snow[3], amount = 0.1)
    ),
    Linewidth_plot = case_when(
      Area == "France" ~ 1.4,
      Area == "Europe" ~ 1.4,
      .default = 0.6
    )
) |> 
  mutate(
    Area = factor(
      Area,
      levels = c(
        sort(setdiff(unique(Area), c("France", "Europe"))),
        "France",
        "Europe"
      )
    )
  ) |> 
  mutate(
    priority = Area %in% c("France", "Europe")
  ) |> 
  arrange(priority)

df_labels <- df_data_FAO %>%
  group_by(Area) %>%
  filter(Year == max(Year)) %>%
  ungroup()

levels(df_data_FAO$Area)

xmax <- max(df_data_FAO$Year)
xmin <- 2000
ymax = 12
ymin = 0

p <- ggplot(
  data = df_data_FAO,
  aes(
    x = Year,
    y = Value,
    color = Color_plot,
    group = Area
  )
)+
  geom_vline(
    xintercept = seq(2000, 2023, by = 1),
    color = "grey91", 
    size = .4
  ) +
  geom_segment(
    data = tibble(y = seq(0, 12, by = 2), x1 = 2000, x2 = 2023),
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    color = "grey91",
    size = .4
  ) +
  geom_line(
    aes(size = Linewidth_plot)
    #linewidth = 1
  )+
  geom_label_repel(
    data = df_labels,
    aes(label = Area_french),
    nudge_x = 4,
    xmin = c(xmax, NA),
    direction = "y",
    hjust = 0,
    fontface = "bold",
    segment.size = .7,
    segment.alpha = .5,
    segment.linetype = "dashed",
    min.segment.length = 0.2,
    box.padding = .4,
    size = 5,
  )+
  scale_x_continuous(
    expand = c(0.01, 0),
    limits = c(xmin, xmax+4),
    breaks = seq(xmin, xmax, by = 2),
    minor_breaks = seq(xmin, xmax, by = 1)
  ) +
  scale_y_continuous(
    breaks = seq(ymin, ymax, by = 2),
    labels = label_comma()
    )+
  scale_color_identity()+
  scale_size_identity()+
  labs(
    y = "Quantité de produits phytopharmaceutiques utilisés\npar surface cultivée (kg/ha)",
    x = "Année"
  )+
  theme_minimal(18)+
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

p

ggsave(
  plot = p, 
  filename = "Ventes_de_pesticides_en_France_et_Europe.png",
  path = path_fig,
  width = 12,
  height = 8
)
