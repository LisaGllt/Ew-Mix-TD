library(here)
source(file = here::here("functions/fun.R"))
f_load_libraries_colors()

# 1. Data importation ----

# Retriving EC50 
df_drc_param <- read.csv(here::here("data/DRC_parameters_cocoons.csv"))
EPX_EC50 <- df_drc_param$coef.drc.4.[3]
IMD_EC50 <- df_drc_param$coef.drc.4.[4]

df_mix_ind <- read_excel(here::here("data/Data_Mix.xlsx"), sheet="ReproSpring2025weight") |> 
  mutate(
    Condition = paste0(Ratio, Line),
    Nb_rep = as.factor(Nb_rep),
    w = as.numeric(w),
    L = w^(1/3),
    TU = Dose_IMD / IMD_EC50 + Dose_EPX/EPX_EC50
  )
df_mix_coc <- read_excel(here::here("data/Data_Mix.xlsx"), sheet="ReproSpring2025coc") |> 
  mutate(
    Condition = paste0(Ratio, Line),
    w_coc = w_coc_tot/(Nb_cocoons-Nb_cocoons_crushed),
    TU = Dose_IMD / IMD_EC50 + Dose_EPX/EPX_EC50
  )

df_conc <- df_mix_coc |> 
  dplyr::select(Condition, Dose_IMD, Dose_EPX)%>%
  distinct(Condition, .keep_all = TRUE)

df_mix_cocsize <- read_excel(here::here("data/Data_Mix.xlsx"), sheet="ReproSpring2025cocsize") |> 
  mutate(
    Condition = paste0(Ratio, Line),
    v_coc = 4/3*pi*Coc_height*Coc_width^2,
  ) |> 
  left_join(df_conc, by = "Condition") |> 
  mutate(TU = Dose_IMD / IMD_EC50 + Dose_EPX/EPX_EC50)

df_mix_cocsize_cosm <- df_mix_cocsize |> 
  aggregate(v_coc ~ Ratio + Line + Nb_rep + TU, FUN = mean) |> 
  mutate(Condition = paste0(Ratio, Line))

df_mix_coc <- df_mix_coc |> 
  left_join(df_mix_cocsize_cosm, by = c("Condition", "Nb_rep", "Ratio", "Line", "TU")) |> 
  mutate(
    Id_Cond = w_coc/v_coc
  )

df_mix_coc_mean <- df_mix_coc |> 
  aggregate(cbind(Nb_cocoons, w_coc, v_coc) ~ Condition + Dose_IMD + Dose_EPX + Ratio + Line + TU, mean)



# Data EC50

df_repro_tot <- read_excel(here::here("data/Data_EC50.xlsx"), sheet="Repro") |> 
  mutate(Dose = as.numeric(Dose)) |> 
  mutate(Dose_f = as.factor(round(Dose, 4))) |> 
  mutate(w=w_tot/Nb_ind) |> 
  mutate(L=w^(1/3))

# Data for EPX
df_repro_EPX <- subset(df_repro_tot, Molec %in% c("Ctrl_2", "EPX")) |> 
  mutate(Molecule = "EPX")
# We only keep earthworms still alive at then end of the experiment
ID_alive_EPX <- subset(df_repro_EPX, t==28 & (!Status=="D"))$ID_cosm
df_repro_EPX_alive <- subset(df_repro_EPX, ID_cosm %in% ID_alive_EPX)

# Data for IMD
df_repro_IMD <- subset(df_repro_tot, Molec %in% c("Ctrl_1", "Ctrl_2", "IMD")) |> 
  mutate(Molecule = "IMD")
# We only keep earthworms still alive at then end of the experiment
ID_alive_IMD <- subset(df_repro_IMD, t==28 & (!Status=="D"))$ID_cosm
df_repro_IMD_alive <- subset(df_repro_IMD, ID_cosm %in% ID_alive_IMD)

df_repro <- rbind(df_repro_EPX, df_repro_IMD)
df_repro_alive <- rbind(df_repro_EPX_alive, df_repro_IMD_alive)
df_repro_alive_mean <- df_repro_alive |> 
  aggregate(Nb_cocoons~Dose+Molecule, mean)

Dose_x <- expand.grid(
  exp(
    seq(
      log(0.000001),log(5000), 
      by=(log(5000)-log(0.000001))/100
    )
  )
)

# Common Ymax and slope
drc.4<- drm(
  Nb_cocoons~Dose, Molecule, 
  data = df_repro_alive,
  type="Poisson",
  fct=LL.4(
    names = c("slope", "Ymin", "Ymax", "EC50"),
    fixed = c(NA, 0, NA, NA)
  ), 
  # Reproduction is expected to reach 0 for large concentrations
  pmodels=data.frame(1, 1, Molecule)
)
# summary(drc.4)

CI.4 <- data.frame(
  Dose = c(Dose_x$Var1, Dose_x$Var1),
  Molecule = c(rep("EPX", length(Dose_x$Var1)), rep("IMD", length(Dose_x$Var1)))
)

pm.4 <- predict(drc.4, newdata=CI.4, interval="confidence")
CI.4$p <- pm.4[,1]
CI.4$pmin <- pm.4[,2]
CI.4$pmax <- pm.4[,3]

Hill.4 <- function(Dose, Molecule){
  
  i <- 0
  if(Molecule == "IMD"){i <- 1}
  
  Y_min <- 0
  slope <- coef(drc.4)[1]
  Y_max <- coef(drc.4)[2]
  EC50 <- coef(drc.4)[3+i]
  
  return(Y_min+(Y_max-Y_min)/(1+exp(slope*(log(Dose)-log(EC50)))))
}
EPX_EC50 <- coef(drc.4)[4]
IMD_EC50 <- coef(drc.4)[4]

x <- seq(0.000001, 5000, 0.01)

df_mix_coc_IMD <- subset(df_mix_coc, Ratio %in% c("I", "N"))
df_mix_coc_EPX <- subset(df_mix_coc, Ratio  %in% c("E", "N"))
df_mix_coc_IMD_mean <- df_mix_coc_IMD |> 
  aggregate(Nb_cocoons ~ Condition+Dose_IMD, mean)
df_mix_coc_EPX_mean <- df_mix_coc_EPX |> 
  aggregate(Nb_cocoons ~ Condition+Dose_EPX, mean)

df_design <- read.csv(here::here("data/Design_mixture.csv"))
df_drc_param <- read.csv(here::here("data/DRC_parameters_cocoons.csv"))

slope <- df_drc_param$coef.drc.4.[1]
Y_max <- df_drc_param$coef.drc.4.[2]
EPX_EC50 <- df_drc_param$coef.drc.4.[3]
IMD_EC50 <- df_drc_param$coef.drc.4.[4]
y_min <- 0

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


grid_C1 <- exp(seq(log(min_dose_EPX), log(max_dose_EPX), by = by_EPX))
grid_C2 <- exp(seq(log(min_dose_IMD), log(max_dose_IMD), by = by_IMD))
grid <- expand.grid(
  x = grid_C1, 
  y = grid_C2
)

# Data set used for the interaction estimation
df_mix_coc_int <- df_mix_coc |> 
  mutate(Response = Nb_cocoons) |> 
  dplyr::select(Dose_EPX, Dose_IMD, Response) |> 
  as.data.frame() |> 
  filter(!(Dose_IMD == 0 & Dose_EPX == 0))

N_obs <- length(df_mix_coc_int$Response)

# 2. Resampling & model fit ----

N_draws <- 10000

# Pick a number btw 1 and N_obs for each draw
min_resampling <- 1
max_resampling <- N_obs

# Draws
draws <- lapply(
  1:N_draws, 
  function(i) {
    n <- sample(min_resampling:max_resampling, 1) 
    lignes <- sample(1:N_obs, size = n, replace = TRUE)  
    list(ID = i, N_resampling = n, Lignes = lignes)
  }
  )

# Convertir en data.frame avec une colonne liste
df_draws <- do.call(rbind, lapply(draws, as.data.frame))
df_draws <- as.data.frame(df_draws)
df_draws$Lignes <- lapply(df_draws$Lignes, unlist) 

str(df_draws)

List_a <- vector("list", N_draws) 
i <- 1

for(draw_i in unique(df_draws$ID)){
  print(i)
  
  Lignes_i <- unlist(subset(df_draws, ID==draw_i)$Lignes)
  
  df_mix_coc_int_i <- df_mix_coc_int |> 
    slice(Lignes_i)
  
  List_a[[i]] <- tryCatch({
    # Code that can fail
    CA_complete_fit_speed(                         
      C_mat=cbind(                                              
        df_mix_coc_int_i$Dose_EPX,                                  
        df_mix_coc_int_i$Dose_IMD                                  
      ),                                                         
      Response=df_mix_coc_int_i$Response,                         
      interact="SA",                                            
      param=param                                                   
    )$a         
  }, error = function(e) {
    NA  # ou NA_real_, NA_character_, etc. selon le type attendu
  })
   i <- i+1
}

quantile(List_a, c(0.025, 0.975))
quantile(List_a, probs = c(0.025, 0.975), na.rm = TRUE)

quantile(unlist(List_a), probs = c(0.025, 0.975), na.rm = TRUE)


