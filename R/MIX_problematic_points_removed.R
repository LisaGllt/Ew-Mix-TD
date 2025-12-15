library(here)
source(file = here::here("functions/fun.R"))
f_load_libraries_colors()

palette <- colorRampPalette(c("#5E81AC", "#88C0D0"))
col_mix <- Nord_polar[4]
pal_ratio <- c("#88C0D0", col_mix, col_mix, col_mix, "#5E81AC",col_mix)
Val_Ctrl <- 1e-4
Label_Ctrl <- "0 (Ctrl)"

# 1. Data ----

# Retriving EC50 
df_drc_param <- read.csv(here::here("data/DRC_parameters_cocoons_really_used.csv"))
EPX_EC50_drc <- df_drc_param$coef.drc.4.[3]
IMD_EC50_drc <- df_drc_param$coef.drc.4.[4]

df_mix_ind <- read_excel(here::here("data/Data_Mix.xlsx"), sheet = "ReproSpring2025weight") |>
  mutate(
    Condition = paste0(Ratio, Line),
    Nb_rep = as.factor(Nb_rep),
    w = as.numeric(w),
    L = w^(1 / 3),
    TU_drc = Dose_IMD / IMD_EC50_drc + Dose_EPX / EPX_EC50_drc,
    Dose_IMD_plot = case_when(
      Dose_IMD == 0 ~ Val_Ctrl,
      !(Dose_IMD == 0) ~ Dose_IMD
    ),
    Dose_EPX_plot = case_when(
      Dose_EPX == 0 ~ Val_Ctrl,
      !(Dose_EPX == 0) ~ Dose_EPX
    ),
    TU_drc_plot = case_when(
      TU_drc == 0 ~ Val_Ctrl,
      !(TU_drc == 0) ~ TU_drc
    )
  )

df_mix_coc <- read_excel(here::here("data/Data_Mix.xlsx"), sheet = "ReproSpring2025coc") |>
  mutate(
    Condition = paste0(Ratio, Line),
    w_coc = w_coc_tot / (Nb_cocoons - Nb_cocoons_crushed),
    TU_drc = Dose_IMD / IMD_EC50_drc + Dose_EPX / EPX_EC50_drc,
    Dose_IMD_plot = case_when(
      Dose_IMD == 0 ~ Val_Ctrl,
      !(Dose_IMD == 0) ~ Dose_IMD
    ),
    Dose_EPX_plot = case_when(
      Dose_EPX == 0 ~ Val_Ctrl,
      !(Dose_EPX == 0) ~ Dose_EPX
    ),
    TU_drc_plot = case_when(
      TU_drc == 0 ~ Val_Ctrl,
      !(TU_drc == 0) ~ TU_drc
    )
  )

df_coc_mean_controls <- subset(df_mix_coc, Dose_IMD == 0 & Dose_EPX == 0) |> 
  aggregate(Nb_cocoons ~ Lot, FUN = mean)

Mean_coc_controls <- mean(subset(df_mix_coc, Dose_IMD == 0 & Dose_EPX == 0)$Nb_cocoons)

#Shift <- 1e-6

df_mix_coc <- df_mix_coc |> 
  mutate(
    r_Nb_cocoons = Nb_cocoons/Mean_coc_controls
  )

df_mix_coc <- df_mix_coc |>
  mutate(
    Trans_r_Nb_cocoons = sqrt(r_Nb_cocoons)
  )

df_conc_mix <- df_mix_coc |> 
  dplyr::select(Condition, Dose_IMD, Dose_EPX)%>%
  distinct(Condition, .keep_all = TRUE)

df_mix_cocsize <- read_excel(here::here("data/Data_Mix.xlsx"), sheet="ReproSpring2025cocsize") |> 
  mutate(
    Condition = paste0(Ratio, Line),
    v_coc = 4/3*pi*(Coc_height/2)*(Coc_width/2)^2,
  ) |> 
  left_join(df_conc_mix, by = "Condition") |> 
  mutate(TU_drc = Dose_IMD / IMD_EC50_drc + Dose_EPX/EPX_EC50_drc)

df_mix_cocsize_cosm <- df_mix_cocsize |> 
  aggregate(v_coc ~ Ratio + Line + Nb_rep + TU_drc, FUN = mean) |> 
  mutate(Condition = paste0(Ratio, Line))

df_mix_coc <- df_mix_coc |> 
  full_join(df_mix_cocsize_cosm, by = c("Condition", "Nb_rep", "Ratio", "Line", "TU_drc")) |> 
  mutate(
    Id_Cond = w_coc/v_coc,
    Condition_f = as.factor(Condition)
  )

df_mix_coc_mean <- df_mix_coc |>
  dplyr::group_by(Condition, Dose_IMD, Dose_EPX, Dose_IMD_plot, Dose_EPX_plot, Ratio, Line, TU_drc) |>
  dplyr::summarise(
    Nb_cocoons = mean(Nb_cocoons, na.rm = TRUE),
    w_coc = mean(w_coc, na.rm = TRUE),
    v_coc = mean(v_coc, na.rm = TRUE),
    r_Nb_cocoons = mean(r_Nb_cocoons, na.rm = TRUE),
    Trans_r_Nb_cocoons = mean(Trans_r_Nb_cocoons, na.rm = TRUE),
    #Nb_cocoons_normalized = mean(Nb_cocoons_normalized, na.rm = TRUE),
    .groups = "drop"
  )

df_mix_coc_IMD_mean <- subset(df_mix_coc_mean, Dose_EPX ==0)
df_mix_coc_EPX_mean <- subset(df_mix_coc_mean, Dose_IMD ==0)


df_mix_coc_sd <- df_mix_coc |>
  dplyr::group_by(Condition, Dose_IMD, Dose_EPX, Dose_IMD_plot, Dose_EPX_plot, Ratio, Line, TU_drc) |>
  dplyr::summarise(
    Nb_cocoons = sd(Nb_cocoons, na.rm = TRUE),
    w_coc = sd(w_coc, na.rm = TRUE),
    v_coc = sd(v_coc, na.rm = TRUE),
    r_Nb_cocoons = sd(r_Nb_cocoons, na.rm = TRUE),
    Trans_r_Nb_cocoons = sd(Trans_r_Nb_cocoons, na.rm = TRUE),
    #Nb_cocoons_normalized = mean(Nb_cocoons_normalized, na.rm = TRUE),
    .groups = "drop"
  )

df_mix_coc_cv <- df_mix_coc_mean |>
  dplyr::rename_with(~ paste0(.x, "_mean"), -c(Condition, Dose_IMD, Dose_EPX, Dose_IMD_plot, Dose_EPX_plot, Ratio, Line, TU_drc)) |>
  dplyr::left_join(
    df_mix_coc_sd |> 
      dplyr::rename_with(~ paste0(.x, "_sd"), -c(Condition, Dose_IMD, Dose_EPX, Dose_IMD_plot, Dose_EPX_plot, Ratio, Line, TU_drc)),
    by = c("Condition", "Dose_IMD", "Dose_EPX", "Dose_IMD_plot", "Dose_EPX_plot", "Ratio", "Line", "TU_drc")
  ) |>
  dplyr::mutate(
    CV_Nb_cocoons = Nb_cocoons_sd / Nb_cocoons_mean,
    CV_w_coc = w_coc_sd / w_coc_mean,
    CV_v_coc = v_coc_sd / v_coc_mean,
    CV_r_Nb_cocoons = r_Nb_cocoons_sd / r_Nb_cocoons_mean,
    CV_Trans_r_Nb_cocoons = Trans_r_Nb_cocoons_sd / Trans_r_Nb_cocoons_mean
  )

df_mix_coc_IMD <- subset(df_mix_coc, Ratio %in% c("I", "N")) |> 
  mutate(
    Dose = Dose_IMD,
    Dose_plot = Dose_IMD_plot,
    Molecule = "IMD"
  ) |> 
  dplyr::select(Dose, Dose_plot, Nb_cocoons, Molecule, Line)

# We do not take the controls into account twice
df_mix_coc_EPX <- subset(df_mix_coc, Ratio  %in% c("E")) |> 
  mutate(
    Dose = Dose_EPX,
    Dose_plot = Dose_EPX_plot,
    Molecule = "EPX"
  ) |>  
  dplyr::select(Dose, Dose_plot, Nb_cocoons, Molecule, Line)

df_mix_coc_EPX_plot <- subset(df_mix_coc, Ratio  %in% c("E", "N")) |>
  mutate(
    Dose = Dose_EPX,
    Dose_plot = Dose_EPX_plot,
    Molecule = "EPX"
  ) |>
  dplyr::select(Dose, Dose_plot, Nb_cocoons, Molecule)

df_mix_coc_single <- rbind(df_mix_coc_IMD, df_mix_coc_EPX)

drc.mix.IE<- drm(
  Nb_cocoons~Dose, Molecule, 
  data = df_mix_coc_single,
  type="Poisson",
  fct=LL.4(
    names = c("slope", "Ymin", "Ymax", "EC50"),
    fixed = c(NA, 0, NA, NA)
  ), 
  # Reproduction is expected to reach 0 for large concentrations
  pmodels=data.frame(1, 1, Molecule)
)

Hill.mix.IE <- function(Dose, Molecule){
  
  i <- ifelse(Molecule == "IMD", 1, 0)
  
  Y_min <- 0
  slope <- coef(drc.mix.IE)[1]
  Y_max <- coef(drc.mix.IE)[2]
  EC50 <- coef(drc.mix.IE)[3+i]
  
  return(Y_min+(Y_max-Y_min)/(1+exp(slope*(log(Dose)-log(EC50)))))
}

df_mix_coc_single <- df_mix_coc_single |> 
  mutate(
    Nb_cocoons_pred = Hill.mix.IE(Dose, Molecule)
  )

LL_drc_Poisson = sum(dpois(df_mix_coc_single$Nb_cocoons, lambda = df_mix_coc_single$Nb_cocoons_pred, log = TRUE))


# 2. Dispersion ----

m.Poisson <- glm(Nb_cocoons ~ Condition_f, family = poisson, data = df_mix_coc)
dispersion <- sum(residuals(m.Poisson, type = "pearson")^2) / m.Poisson$df.residual
dispersion
library(AER)
dispersiontest(m.Poisson)


#m.QPoisson <- survey::svyglm(Nb_cocoons ~ Condition_f, family = quasipoisson, data = df_mix_coc)
m.NB <- glm.nb(Nb_cocoons ~ Condition_f, data = df_mix_coc)

# 1. Calculer moyenne et variance par groupe
library(dplyr)

variance_mean <- df_mix_coc %>%
  group_by(Condition_f) %>%
  summarise(
    mean = mean(Nb_cocoons),
    variance = var(Nb_cocoons),
    n = n()
  )

theta_nb <- summary(m.NB)$theta

# 3. Créer les données pour les courbes théoriques
mean_range <- seq(0, max(variance_mean$mean) * 1.1, length.out = 100)

theoretical <- data.frame(
  mean = rep(mean_range, 3),
  variance = c(
    mean_range,  # Poisson
    dispersion * mean_range,  # Quasi-Poisson
    mean_range + mean_range^2 / theta_nb  # Negative Binomial
  ),
  model = rep(c("Poisson", "Quasi-Poisson", "Binomiale Negative"), 
              each = length(mean_range))
)

# 4. Créer le graphique
ggplot() +
  # Lignes théoriques
  geom_line(data = theoretical, 
            aes(x = mean, y = variance, color = model, linetype = model), 
            linewidth = 1.2) +
  # Points observés
  geom_point(data = variance_mean, 
             aes(x = mean, y = variance), 
             size = 4, color = "black", alpha = 0.5) +
  # Labels des conditions
  geom_text(data = variance_mean, 
            aes(x = mean, y = variance, label = Condition_f), 
            vjust = -1.2, hjust = 0.5, size = 3.5, fontface = "bold") +
  # Mise en forme
  labs(x = "Moyenne (μ)", 
       y = "Variance (σ²)", 
       title = "Relation Variance-Moyenne : Choix du modèle de distribution",
       subtitle = "Les points noirs représentent les valeurs observées par condition",
       color = "Distribution théorique",
       linetype = "Distribution théorique") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30")
  ) +
  scale_color_manual(values = c(
    "Poisson" = "gray50", 
    "Quasi-Poisson" = "#E41A1C", 
    "Binomiale Negative" = "#377EB8"
  )) +
  scale_linetype_manual(values = c(
    "Poisson" = "dashed",
    "Quasi-Poisson" = "solid",
    "Binomiale Negative" = "solid"
  )) +
  guides(color = guide_legend(nrow = 3), 
         linetype = guide_legend(nrow = 3))

# Afficher aussi le tableau des valeurs
print(variance_mean)

# 3. Préparation des entrées des modèles MIX ----

# Points problématiques enlevés !!!!!!!!!!!!!!!!!!!!
df_mix_coc_filter <- df_mix_coc |> 
  filter(!Condition %in% c("F2", "I2", "H1", "E5"))

df_mix_coc_IMD_filter <- subset(df_mix_coc_filter, Ratio %in% c("I", "N")) |> 
  mutate(
    Dose = Dose_IMD,
    Dose_plot = Dose_IMD_plot,
    Molecule = "IMD"
  ) |> 
  dplyr::select(Dose, Dose_plot, Nb_cocoons, Molecule, Line)

# We do not take the controls into account twice
df_mix_coc_EPX_filter <- subset(df_mix_coc_filter, Ratio  %in% c("E")) |> 
  mutate(
    Dose = Dose_EPX,
    Dose_plot = Dose_EPX_plot,
    Molecule = "EPX"
  ) |>  
  dplyr::select(Dose, Dose_plot, Nb_cocoons, Molecule, Line)

df_mix_coc_EPX_filter_plot <- subset(df_mix_coc_filter, Ratio  %in% c("E", "N")) |>
  mutate(
    Dose = Dose_EPX,
    Dose_plot = Dose_EPX_plot,
    Molecule = "EPX"
  ) |>
  dplyr::select(Dose, Dose_plot, Nb_cocoons, Molecule)

df_mix_coc_single_filter <- rbind(df_mix_coc_IMD_filter, df_mix_coc_EPX_filter) 

df_design <- read.csv(here::here("data/Design_mixture.csv"))

drc.mix.IE.filter<- drm(
  Nb_cocoons~Dose, Molecule, 
  data = df_mix_coc_single_filter,
  type="Poisson",
  fct=LL.4(
    names = c("slope", "Ymin", "Ymax", "EC50"),
    fixed = c(NA, 0, NA, NA)
  ), 
  # Reproduction is expected to reach 0 for large concentrations
  pmodels=data.frame(1, 1, Molecule)
)

Hill.mix.IE.filter <- function(Dose, Molecule){
  
  i <- ifelse(Molecule == "IMD", 1, 0)
  
  Y_min <- 0
  slope <- coef(drc.mix.IE.filter)[1]
  Y_max <- coef(drc.mix.IE.filter)[2]
  EC50 <- coef(drc.mix.IE.filter)[3+i]
  
  return(Y_min+(Y_max-Y_min)/(1+exp(slope*(log(Dose)-log(EC50)))))
}

df_mix_coc_single_filter <- df_mix_coc_single_filter |> 
  mutate(
    Nb_cocoons_pred = Hill.mix.IE.filter(Dose, Molecule)
  )

LL_drc_Poisson_filter = sum(dpois(df_mix_coc_single_filter$Nb_cocoons, lambda = df_mix_coc_single_filter$Nb_cocoons_pred, log = TRUE))


df_mix_coc_int <- df_mix_coc |> 
  mutate(Response = Nb_cocoons) |> 
  dplyr::select(Dose_EPX, Dose_IMD, Response, Line, Ratio, TU_drc, Nb_rep) |> 
  as.data.frame() |> 
  filter(!(Dose_IMD == 0 & Dose_EPX == 0))

df_mix_coc_int_single <- df_mix_coc_int |> 
  filter(Dose_EPX == 0 | Dose_IMD == 0)

df_mix_coc_int_filter <- df_mix_coc_filter |> 
  mutate(Response = Nb_cocoons) |> 
  dplyr::select(Dose_EPX, Dose_IMD, Response, Line, Ratio, TU_drc, Nb_rep) |> 
  as.data.frame() |> 
  filter(!(Dose_IMD == 0 & Dose_EPX == 0))

df_mix_coc_int_single_filter <- df_mix_coc_int_filter |> 
  filter(Dose_EPX == 0 | Dose_IMD == 0)


C_mat = cbind(df_design$EPX,df_design$IMD)
Max = coef(drc.mix.IE)[2]
Slopes = c(coef(drc.mix.IE)[1], coef(drc.mix.IE)[1])
Ec50s = c(coef(drc.mix.IE)[3], coef(drc.mix.IE)[4])

Max_filter = coef(drc.mix.IE)[2]
Slopes_filter = c(coef(drc.mix.IE)[1], coef(drc.mix.IE)[1])
Ec50s_filter = c(coef(drc.mix.IE)[3], coef(drc.mix.IE)[4])

param <- data.frame(Slopes=Slopes, Max=Max, Ec50s=Ec50s)
param_filter <- data.frame(Slopes=Slopes_filter, Max=Max_filter, Ec50s=Ec50s_filter)

min_dose_EPX <- min(subset(df_design, !EPX==0)$EPX)
max_dose_EPX <- max(subset(df_design, !EPX==0)$EPX)

min_dose_IMD <- min(subset(df_design, !IMD==0)$IMD)
max_dose_IMD <- max(subset(df_design, !IMD==0)$IMD)

by_EPX <- (log(max_dose_EPX)-log(min_dose_EPX))/100
by_IMD <- (log(max_dose_IMD)-log(min_dose_IMD))/100

# Créer une grille pour x et y en utilisant expand.grid
grid_C1 <- exp(seq(log(min_dose_EPX), log(max_dose_EPX), by = by_EPX))
grid_C2 <- exp(seq(log(min_dose_IMD), log(max_dose_IMD), by = by_IMD))
grid <- expand.grid(
  x = grid_C1, 
  y = grid_C2
)

# 4. Fit des modèles ----

signif_LL <- 5

########### CA #######################################

CA_fit <- data.frame(Error = NA)
CA_fit$Error <- CA_complete2_Poisson(
  C_mat=cbind(                                                    # <1>
    df_mix_coc_int$Dose_EPX,                                      # <1>
    df_mix_coc_int$Dose_IMD                                       # <1>
  ),                                                            # <1>
  Response=df_mix_coc_int$Response,                               # <1>
  Max = Max,
  Slopes = Slopes,
  Ec50s = Ec50s,
  interact="none",                                                # <1>
  multicore = TRUE, 
  mc.cores = 4
)

CA_fit$Type <- "CA"
CA_fit$Ref <- "CA"

LogLikelihood_CA <- -CA_fit$Error                              
CA_fit$a <- NA
CA_fit$b <- NA
# Printing results
LogLikelihood_CA_print <- signif(LogLikelihood_CA, signif_LL)


CA_fit_filter <- data.frame(Error = NA)
CA_fit_filter$Error <- CA_complete2_Poisson(
  C_mat=cbind(                                                    # <1>
    df_mix_coc_int$Dose_EPX,                                      # <1>
    df_mix_coc_int$Dose_IMD                                       # <1>
  ),                                                            # <1>
  Response=df_mix_coc_int$Response,                               # <1>
  Max = Max_filter,
  Slopes = Slopes_filter,
  Ec50s = Ec50s_filter,
  interact="none",                                                # <1>
  multicore = TRUE, 
  mc.cores = 4
)

CA_fit_filter$Type <- "CA"
CA_fit_filter$Ref <- "CA"

LogLikelihood_CA_filter <- -CA_fit_filter$Error                              
CA_fit_filter$a <- NA
CA_fit_filter$b <- NA
# Printing results
LogLikelihood_CA_filter_print <- signif(LogLikelihood_CA_filter, signif_LL)

################# IA ####################################


IA_fit <- data.frame(Error = NA)
IA_fit$Error <- IA_complete2_Poisson(
  C_mat=cbind(                                                    # <1>
    df_mix_coc_int$Dose_EPX,                                      # <1>
    df_mix_coc_int$Dose_IMD                                       # <1>
  ),                                                            # <1>
  Response=df_mix_coc_int$Response,                               # <1>
  Max = Max,
  Slopes = Slopes,
  Ec50s = Ec50s,
  interact="none"
)

IA_fit$Type <- "IA"
IA_fit$Ref <- "IA"

LogLikelihood_IA <- -IA_fit$Error                              
IA_fit$a <- NA
IA_fit$b <- NA
# Printing results
LogLikelihood_IA_print <- signif(LogLikelihood_IA, signif_LL)


IA_fit_filter <- data.frame(Error = NA)
IA_fit_filter$Error <- IA_complete2_Poisson(
  C_mat=cbind(                                                    # <1>
    df_mix_coc_int$Dose_EPX,                                      # <1>
    df_mix_coc_int$Dose_IMD                                       # <1>
  ),                                                            # <1>
  Response=df_mix_coc_int$Response,                               # <1>
  Max = Max_filter,
  Slopes = Slopes_filter,
  Ec50s = Ec50s_filter,
  interact="none"
)

IA_fit_filter$Type <- "IA"
IA_fit_filter$Ref <- "IA"

LogLikelihood_IA_filter <- -IA_fit_filter$Error                              
IA_fit_filter$a <- NA
IA_fit_filter$b <- NA
# Printing results
LogLikelihood_IA_filter_print <- signif(LogLikelihood_IA_filter, signif_LL)

######## CASA ##################################################################

CASA_fit<-CA_complete_fit_speed(                                    # <1>
  C_mat=cbind(                                                    # <1>
    df_mix_coc_int$Dose_EPX,                                      # <1>
    df_mix_coc_int$Dose_IMD                                       # <1>
  ),                                                            # <1>
  Response=df_mix_coc_int$Response,                               # <1>
  interact="SA",                                                  # <1>
  param=param,                                               # <1>
  error_type = "Poisson",
  multicore = TRUE,                                               # <1>
  mc.cores = 4                                                    # <1>
) |>                                                                # <1>
  as.data.frame()
CASA_fit$Type <- "SA"
CASA_fit$Ref <- "CA"

LogLikelihood_CASA <- -CASA_fit$Error                              
a_CASA <- CASA_fit$a
CASA_fit$b <- NA
# Printing results
a_CASA_print <- signif(a_CASA, 3)
LogLikelihood_CASA_print <- signif(LogLikelihood_CASA, signif_LL)

CASA_fit_filter<-CA_complete_fit_speed(                                    # <1>
  C_mat=cbind(                                                    # <1>
    df_mix_coc_int$Dose_EPX,                                      # <1>
    df_mix_coc_int$Dose_IMD                                       # <1>
  ),                                                            # <1>
  Response=df_mix_coc_int$Response,                               # <1>
  interact="SA",                                                  # <1>
  param=param_filter,                                                   # <1>
  error_type = "Poisson",
  multicore = TRUE,                                               # <1>
  mc.cores = 4                                                    # <1>
) |>                                                                # <1>
  as.data.frame()
CASA_fit_filter$Type <- "SA"
CASA_fit_filter$Ref <- "CA"

LogLikelihood_CASA_filter <- -CASA_fit_filter$Error                              
a_CASA_filter <- CASA_fit_filter$a
CASA_fit_filter$b <- NA
# Printing results
a_CASA_filter_print <- signif(a_CASA_filter, 3)
LogLikelihood_CASA_filter_print <- signif(LogLikelihood_CASA_filter, signif_LL)

######## IASA ##################################################################

IASA_fit<-IA_complete_fit_speed(                                    # <1>
  C_mat=cbind(                                                    # <1>
    df_mix_coc_int$Dose_EPX,                                      # <1>
    df_mix_coc_int$Dose_IMD                                       # <1>
  ),                                                            # <1>
  Response=df_mix_coc_int$Response,                               # <1>
  interact="SA",                                                  # <1>
  param=param,                                                 # <1>
  error_type = "Poisson",
  multicore = TRUE,                                               # <1>
  mc.cores = 4                                                    # <1>
) |>                                                                # <1>
  as.data.frame()
IASA_fit$Type <- "SA"
IASA_fit$Ref <- "IA"

LogLikelihood_IASA <- -IASA_fit$Error                              
a_IASA <- IASA_fit$a
IASA_fit$b <- NA
# Printing results
a_IASA_print <- signif(a_IASA, 3)
LogLikelihood_IASA_print <- signif(LogLikelihood_IASA, signif_LL)

IASA_fit_filter<-IA_complete_fit_speed(                                    # <1>
  C_mat=cbind(                                                    # <1>
    df_mix_coc_int$Dose_EPX,                                      # <1>
    df_mix_coc_int$Dose_IMD                                       # <1>
  ),                                                            # <1>
  Response=df_mix_coc_int$Response,                               # <1>
  interact="SA",                                                  # <1>
  param=param_filter,                                                   # <1>
  error_type = "Poisson",
  multicore = TRUE,                                               # <1>
  mc.cores = 4                                                    # <1>
) |>                                                                # <1>
  as.data.frame()
IASA_fit_filter$Type <- "SA"
IASA_fit_filter$Ref <- "IA"

LogLikelihood_IASA_filter <- -IASA_fit_filter$Error                              
a_IASA_filter <- IASA_fit_filter$a
IASA_fit_filter$b <- NA
# Printing results
a_IASA_filter_print <- signif(a_IASA_filter, 3)
LogLikelihood_IASA_filter_print <- signif(LogLikelihood_IASA_filter, signif_LL)





