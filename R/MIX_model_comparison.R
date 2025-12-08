library(tidyverse)
library(readxl)
library(drc)
library(here)
library(parallel)
#library(gl)

source(file = here::here("functions/fun.R"))

# 1. Importation des données ----

# Mixture data
Val_Ctrl <- 1e-4
Label_Ctrl <- "0 (Ctrl)"

df_mix_ind <- read_excel(here::here("data/Data_Mix.xlsx"), sheet = "ReproSpring2025weight") |>
  mutate(
    Condition = paste0(Ratio, Line),
    Nb_rep = as.factor(Nb_rep),
    w = as.numeric(w),
    L = w^(1 / 3),
    Dose_IMD_plot = case_when(
      Dose_IMD == 0 ~ Val_Ctrl,
      !(Dose_IMD == 0) ~ Dose_IMD
    ),
    Dose_EPX_plot = case_when(
      Dose_EPX == 0 ~ Val_Ctrl,
      !(Dose_EPX == 0) ~ Dose_EPX
    )
  )

df_mix_coc <- read_excel(here::here("data/Data_Mix.xlsx"), sheet = "ReproSpring2025coc") |>
  mutate(
    Condition = paste0(Ratio, Line),
    w_coc = w_coc_tot / (Nb_cocoons - Nb_cocoons_crushed),
    Dose_IMD_plot = case_when(
      Dose_IMD == 0 ~ Val_Ctrl,
      !(Dose_IMD == 0) ~ Dose_IMD
    ),
    Dose_EPX_plot = case_when(
      Dose_EPX == 0 ~ Val_Ctrl,
      !(Dose_EPX == 0) ~ Dose_EPX
    )
  )

df_coc_mean_controls <- subset(df_mix_coc, Dose_IMD == 0 & Dose_EPX == 0) |> 
  aggregate(Nb_cocoons ~ Lot, FUN = mean)

Mean_coc_controls <- mean(subset(df_mix_coc, Dose_IMD == 0 & Dose_EPX == 0)$Nb_cocoons)

#Shift <- 1e-6

df_mix_coc_mean <- df_mix_coc |>
  dplyr::group_by(Condition, Dose_IMD, Dose_EPX, Dose_IMD_plot, Dose_EPX_plot, Ratio, Line) |>
  dplyr::summarise(
    Nb_cocoons = mean(Nb_cocoons, na.rm = TRUE),
    w_coc = mean(w_coc, na.rm = TRUE),
    .groups = "drop"
  )

df_mix_coc_IMD_mean <- subset(df_mix_coc_mean, Dose_EPX ==0)
df_mix_coc_EPX_mean <- subset(df_mix_coc_mean, Dose_IMD ==0)


# 2. Etablissement des dose-réponse ----

df_mix_coc_IMD <- subset(df_mix_coc, Ratio %in% c("I", "N")) |> 
  mutate(
    Dose = Dose_IMD,
    Dose_plot = Dose_IMD_plot,
    Molecule = "IMD"
  ) |> 
  dplyr::select(Dose, Dose_plot, Nb_cocoons, Molecule)

# We do not take the controls into account twice
df_mix_coc_EPX <- subset(df_mix_coc, Ratio  %in% c("E")) |> 
  mutate(
    Dose = Dose_EPX,
    Dose_plot = Dose_EPX_plot,
    Molecule = "EPX"
  ) |>  
  dplyr::select(Dose, Dose_plot, Nb_cocoons, Molecule)

df_mix_coc_EPX_plot <- subset(df_mix_coc, Ratio  %in% c("E", "N")) |> 
  mutate(
    Dose = Dose_EPX,
    Dose_plot = Dose_EPX_plot,
    Molecule = "EPX"
  ) |>  
  dplyr::select(Dose, Dose_plot, Nb_cocoons, Molecule)

df_mix_coc_single <- rbind(df_mix_coc_IMD, df_mix_coc_EPX)

df_mix_coc_single_mean <- df_mix_coc_single |> 
  aggregate(Nb_cocoons ~ Dose+Molecule, mean)

Dose_x <- expand.grid(
  exp(
    seq(
      log(Val_Ctrl),log(5000),
      by=(log(5000)-log(Val_Ctrl))/100
    )
  )
)

# Common Ymax and slope
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

CI.mix.IE <- data.frame(
  Dose = c(Dose_x$Var1, Dose_x$Var1),
  Molecule = c(rep("EPX", length(Dose_x$Var1)), rep("IMD", length(Dose_x$Var1)))
)

pm.mix.IE <- predict(drc.mix.IE, newdata=CI.mix.IE, interval="confidence")
CI.mix.IE$p <- pm.mix.IE[,1]
CI.mix.IE$pmin <- pm.mix.IE[,2]
CI.mix.IE$pmax <- pm.mix.IE[,3]

Hill.mix.IE <- function(Dose, Molecule){
  
  i <- 0
  if(Molecule == "IMD"){i <- 1}
  
  Y_min <- 0
  slope <- coef(drc.mix.IE)[1]
  Y_max <- coef(drc.mix.IE)[2]
  EC50 <- coef(drc.mix.IE)[3+i]
  
  return(Y_min+(Y_max-Y_min)/(1+exp(slope*(log(Dose)-log(EC50)))))
}
EPX_EC50_mix.IE <- coef(drc.mix.IE)[3]
IMD_EC50_mix.IE <- coef(drc.mix.IE)[4]

df_mix_coc <- df_mix_coc |> 
  mutate(
    TU_mix = Dose_IMD/IMD_EC50_mix.IE + Dose_EPX/EPX_EC50_mix.IE,
    TU_mix_plot = case_when(
      TU_mix == 0 ~ Val_Ctrl,
      !(TU_mix == 0) ~ TU_mix
    )
  )

df_mix_coc_mean <- df_mix_coc_mean |> 
  mutate(
    TU_mix = Dose_IMD/IMD_EC50_mix.IE + Dose_EPX/EPX_EC50_mix.IE,
    TU_mix_plot = case_when(
      TU_mix == 0 ~ Val_Ctrl,
      !(TU_mix == 0) ~ TU_mix
    )
  )

summary(drc.mix.IE)

est <- coef(drc.mix.IE)
ci <- confint(drc.mix.IE)
params_tab <- data.frame(
  estimate  = unname(est),
  CI_lower  = ci[,1],
  CI_upper  = ci[,2]
)

params_tab

# 3. Fit des modèles ----

df_design <- read.csv(here::here("data/Design_mixture.csv"))

slope <- coef(drc.mix.IE)[1]
Y_max <- coef(drc.mix.IE)[2]
EPX_EC50 <- coef(drc.mix.IE)[3]
IMD_EC50 <- coef(drc.mix.IE)[4]
Y_min <- 0

C_mat = cbind(df_design$EPX,df_design$IMD)
Max = Y_max
Slopes = c(slope, slope)
Ec50s = c(EPX_EC50, IMD_EC50)

param <- data.frame(Slopes=Slopes, Max=Max, Ec50s=Ec50s)

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

# Données utilisée pour le fit
df_mix_coc_int <- df_mix_coc |> 
  mutate(Response = Nb_cocoons) |> 
  dplyr::select(Dose_EPX, Dose_IMD, Response, Line, Ratio, TU_mix) |> 
  as.data.frame() |> 
  filter(!(Dose_IMD == 0 & Dose_EPX == 0))

## 3.1 Modèle CA ----

signif_LL <- 5 # Significative digits for LL printing

######## CA ##################################################################

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

LogLikelihood_CA <- CA_fit$Error 

CA_fit$a <- NA
CA_fit$b <- NA
CA_fit$Type <- "CA"
CA_fit$Ref <- "CA"

CA_fit <- CA_fit |> 
  as.data.frame()

## 3.2 Modèle IA ----

######## IA ##################################################################

IA_fit <- data.frame(Error = NA)
IA_fit$Error <- IA_complete2_Poisson(
  C_mat=cbind(                                                  
    df_mix_coc_int$Dose_EPX,                                    
    df_mix_coc_int$Dose_IMD                                      
  ),                                                          
  Response=df_mix_coc_int$Response,                              
  Max = Max, 
  Slopes = Slopes,
  Ec50s = Ec50s,
  interact="none"
)

LogLikelihood_IA <- IA_fit$Error 

IA_fit$a <- NA
IA_fit$b <- NA
IA_fit$Type <- "IA"
IA_fit$Ref <- "IA"

IA_fit <- IA_fit |> 
  as.data.frame()

## 3.3 Anova ----

df_mix_coc_int_anova <- df_mix_coc_int |> 
  mutate(
    Dose_EPX_f = as.factor(Dose_EPX),
    Dose_IMD_f = as.factor(Dose_IMD)
  )
# The anova model 
glm_anova <- glm(Response ~ Dose_EPX_f + Dose_IMD_f, family = poisson(link = "log"), data = df_mix_coc_int_anova)
LogLikelihood_anova <- -as.numeric(logLik(glm_anova)) # Log-vraissemblance négative
k_anova <- attr(logLik(glm_anova), "df")   # nombre de paramètres (inclut intercept)

# 4. Comparaisons ----

## 4.1 Anova vs. CA ----

LL1 <- LogLikelihood_anova    # log-vraisemblance négative du modèle CA
LL0 <- LogLikelihood_CA    # log-vraisemblance négative du modèle SA (CA reference)

p1 <- k_anova         
p0 <- 4     

D <- -2 * (LL1 - LL0)
p_val <- pchisq(D, df = p1 - p0, lower.tail = FALSE)

cat("Comparison Anova/CA - D =", D, "p =", p_val, "\n")

## 4.2 Anova vs IA ----

LL1 <- LogLikelihood_anova    # log-vraisemblance négative du modèle CA
LL0 <- LogLikelihood_IA    # log-vraisemblance négative du modèle SA (CA reference)

p1 <- k_anova         
p0 <- 4     

D <- -2 * (LL1 - LL0)
p_val <- pchisq(D, df = p1 - p0, lower.tail = FALSE)

cat("Comparison Anova/IA - D =", D, "p =", p_val, "\n")



