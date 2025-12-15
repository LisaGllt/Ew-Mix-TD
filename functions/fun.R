# General setup ----

f_load_libraries <- function() {
  # 📦 Data Manipulation
  library(tidyverse) # Collection de packages pour manipulation et visualisation des données
  library(here) # Gestion des chemins de fichiers
  library(readxl) # Lecture des fichiers Excel
  library(reshape2) # Restructuration des données
  library(DT) # Génération de tables interactives
  library(knitr) # Génération de rapports dynamiques en RMarkdown/Quarto

  # 🎨 Visualization
  library(ggplot2) # Visualisation de données
  library(ggthemes) # Thèmes pour ggplot2
  library(ggdist) # Distribution et incertitude
  library(ggsci) # Palettes de couleurs scientifiques
  library(viridis) # Palette de couleurs perceptuellement uniforme
  library(wesanderson) # Palette de couleurs artistiques
  library(RColorBrewer) # Palettes de couleurs prédéfinies
  library(nord) # Palettes de couleurs inspirées du Nord
  library(plotly) # Graphiques interactifs
  library(ggiraph) # Graphiques interactifs pour ggplot2
  library(ggrepel) # Étiquettes non chevauchantes sur ggplot2
  library(patchwork) # Combinaison de plusieurs ggplots
  library(gridExtra) # Arrangements de graphiques en grille
  library(grid) # Outils de mise en page graphique
  library(ggbreak) # Briser les axes dans ggplot2
  library(ggtext) # Formatage avancé de texte dans ggplot2
  library(kableExtra) # Mise en forme avancée des tables
  library(flextable)
  library(gt) # Tables rendering
  library(processx)
  library(metR)
  library(rnaturalearth)
  library(sf)
  library(scales)
  library(colorspace)

  # 📊 Statistical Modeling & Bayesian Analysis
  library(brms) # Modélisation bayésienne avec Stan
  library(rstan) # Interface R pour Stan
  library(cmdstanr) # Interface alternative pour Stan (CmdStan)
  library(tidybayes) # Manipulation et visualisation des résultats bayésiens
  library(ggmcmc) # Diagnostics des chaînes MCMC
  library(rethinking) # Modélisation bayésienne avancée
  library(priorsense) # Analyse de sensibilité des priors
  library(coda)
  library(ggmcmc)
  library(bayesnec)
  library(truncnorm)

  # 🔬 Regression & Hypothesis Testing
  library(car) # Tests statistiques et régressions avancées
  library(nlstools) # Outils pour modèles non linéaires
  library(lsmeans) # Comparaisons post-hoc
  library(ggpubr) # Outils pour publications scientifiques
  library(marginaleffects) # Effets marginaux des modèles
  library(brglm2) # Régressions logistiques biais-réduits
  library(multcomp)

  # ⚙️ Computational Tools & Parallelization
  library(parallel) # Calcul parallèle
  library(deSolve) # Équations différentielles
  library(tmvtnorm) # Distribution normale tronquée multivariée
  library(fdrtool) # Faux taux de découverte (FDR)
  library(drc) # Modélisation de réponses aux doses

  # 🛠️ Model Evaluation & Performance
  library(easystats) # Outils pour statistiques et modèles
  library(performance) # Diagnostics et évaluation de modèles
  library(modelsummary) # Résumé des modèles statistiques
  library(plan) # Planification de l'exécution des tâches

  # 🖋️ Math & LaTeX Support
  library(latex2exp) # Expressions LaTeX dans ggplot2
  library(extrafont) # Gestion des polices pour ggplot2

  library(data.table)
}


f_load_colors <- function() {
  col_blue <<- "#5E81AC"
  col_red <<- "#f42404"

  pal_blue <- c("#5E81AC", "#7F9DC4", "#A0C1D9", "#DCE9F2")
  pal_red <<- c("#f42404", "#F65E4B", "#F6876D", "#FBD3D0")

  Nord_frost <<- nord(palette = "frost")
  Nord_aurora <<- nord(palette = "aurora")
  Nord_polar <<- nord(palette = "polarnight")
  Nord_snow <<- nord(palette = "snowstorm")

  pal_col <<- c(Nord_aurora[1], Nord_frost[4])
  col_EPX <<- Nord_frost[2]
  col_IMD <<- Nord_frost[4]
  col_elim <<- Nord_aurora[4]
  col_uptake <<- Nord_aurora[1]
  shape_IMD <<- 16
  shape_EPX <<- 15
  col_Molecule <<- c(col_EPX, col_IMD)
  col_molec <<- c(col_EPX, col_IMD)
  sizetitle <<- 12

  shape_Molecule <<- c(shape_EPX, shape_IMD)
  shape_molec <<- c(shape_EPX, shape_IMD)

  set.seed(121212)
}

f_load_libraries_colors <- function() {
  f_load_libraries()
  f_load_colors()
}

# Load experimental data ----

f_read_data_expA_repro <- function(){
  
  # Value given to controls for plots
  Val_Ctrl <- 1e-4
  Label_Ctrl <- "0 (Ctrl)"
  
  df_repro_tot <- read_excel(
    here::here("data/Data_EC50.xlsx"), 
    sheet="Repro"
  ) |> 
    mutate(                                              
      Dose = as.numeric(Dose),                            
      Dose_f = as.factor(round(Dose, 4)),           
      Dose_plot = case_when(                          
        Dose == 0 ~ Val_Ctrl,                          
        !(Dose == 0) ~ Dose                              
      ),                                              
      w=w_tot/Nb_ind,                                  
      L=w^(1/3)                                     
    )                                                
  
  df_controls <- subset(df_repro_tot, Molec %in% c("Ctrl_1", "Ctrl_2")) |> 
    aggregate(Nb_cocoons ~ Molec, FUN = mean)
  
  df_repro_EPX <- subset(df_repro_tot, Molec %in% c("Ctrl_2", "EPX")) |>   
    mutate(
      Molecule = "EPX",
      r_Nb_cocoons = Nb_cocoons/subset(df_controls,Molec == "Ctrl_2")$Nb_cocoons*100
    )                                                      
  ID_alive_EPX <- subset(df_repro_EPX, t==28 & (!Status=="D"))$ID_cosm    
  df_repro_EPX_alive <- subset(df_repro_EPX, ID_cosm %in% ID_alive_EPX)   
  
  df_repro_IMD <- subset(df_repro_tot, Molec %in% c("Ctrl_1", "Ctrl_2", "IMD")) |> 
    mutate(
      Molecule = "IMD",
      r_Nb_cocoons = case_when(
        Lot == "A" ~ Nb_cocoons/subset(df_controls,Molec == "Ctrl_1")$Nb_cocoons,
        Lot == "B" ~ Nb_cocoons/subset(df_controls,Molec == "Ctrl_2")$Nb_cocoons
      )*100
    )                                                      
  
  ID_alive_IMD <- subset(df_repro_IMD, t==28 & (!Status=="D"))$ID_cosm            
  df_repro_IMD_alive <- subset(df_repro_IMD, ID_cosm %in% ID_alive_IMD)           
  
  df_repro <- rbind(df_repro_EPX, df_repro_IMD)                                   
  df_repro_alive <- rbind(df_repro_EPX_alive, df_repro_IMD_alive)                 
  
  df_repro_alive_mean <- df_repro_alive |>
    group_by(Dose, Molecule, Dose_plot) |> 
    summarise(
      r_Nb_cocoons = mean(r_Nb_cocoons, na.rm = TRUE),
      Nb_cocoons = mean(Nb_cocoons, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_res <- list(
    df_expA_repro = df_repro,
    df_expA_repro_alive = df_repro_alive,
    df_expA_repro_alive_mean = df_repro_alive_mean
    
  )
  
  return(df_res)
  
}

f_read_data_expA_hatchlings <- function(){
  
  # Value given to controls for plots
  Val_Ctrl <- 1e-4
  Label_Ctrl <- "0 (Ctrl)"
  
  df_repro_tot <- read_excel(
    here::here("data/Data_EC50.xlsx"), 
    sheet="Repro"
  ) |> 
    mutate(                                              
      Dose = as.numeric(Dose),                            
      Dose_f = as.factor(round(Dose, 4)),           
      Dose_plot = case_when(                          
        Dose == 0 ~ Val_Ctrl,                          
        !(Dose == 0) ~ Dose                              
      ),                                              
      w=w_tot/Nb_ind,                                  
      L=w^(1/3)                                     
    )  
  
  df_hatchlings_tot <- read_excel(here::here("data/Data_EC50.xlsx"), sheet="Hatchlings") |>
    mutate(ID_cosm = ID_cosm) |> 
    group_by(ID_cosm) %>%
    arrange(Date) %>%  # S'assurer que les dates sont triées
    mutate(Nb_juv_cum = cumsum(Nb_juv))  # Calcul du cumul
  
  t0_dates <- df_repro_tot %>%
    filter(t == 0) %>%
    dplyr::select(ID_cosm, Date) %>%
    rename(t0_date = Date)
  
  nb_cocoons_t28 <- df_repro_tot %>%
    filter(t == 28) %>%
    dplyr::select(ID_cosm, Nb_cocoons)
  
  
  df_hatchlings_tot <- df_hatchlings_tot |> 
    left_join(t0_dates, by = "ID_cosm") |> 
    mutate(t = as.numeric(difftime(Date, t0_date, units = "days"))) |>   # Calcul du temps relatif
    left_join(nb_cocoons_t28, by = "ID_cosm") |> 
    mutate(Percent_Juv = (Nb_juv_cum / Nb_cocoons) * 100) |>   # Calcul du pourcentage
    left_join(
      df_expA_repro_tot |>  
        dplyr::select(ID_cosm, Dose) |>  
        distinct(), 
      by = "ID_cosm"
    )
  
  df_hatchlings_f <- df_hatchlings_tot |> 
    group_by(ID_cosm) |>           # Grouper par individu
    slice_max(Nb_juv_cum, n = 1) |>   # Sélectionner la ligne avec le maximum
    ungroup() |>                       # Supprimer le groupement
    aggregate(Nb_juv_cum ~ Dose+Molec, FUN = mean) |> 
    mutate(Molecule = Molec)
  
  return(df_hatchlings_f)
}

f_read_data_expA_growth <- function(){
  # Value given to controls for plots
  Val_Ctrl <- 1e-4
  Label_Ctrl <- "0 (Ctrl)"
  
  df_growth_tot <- read_excel(
    here::here("data/Data_EC50.xlsx"), 
    sheet="Growth"
  ) |>                                          
    mutate(                                       
      Dose = as.numeric(Dose),                    
      Dose_f = as.factor(round(Dose, 4)),         
      Dose_plot = case_when(                      
        Dose == 0 ~ Val_Ctrl,                     
        !(Dose == 0) ~ Dose                       
      ),                                          
      L=w^(1/3)                                   
    )                                          

  
  df_growth_EPX <- subset(df_growth_tot, Molec %in% c("Ctrl_2", "EPX")) |> 
    mutate(Molecule = "EPX")                                          
  ID_alive_EPX <- subset(df_growth_EPX, t==28 & (!Status=="D"))$ID      
  df_growth_EPX_alive <- subset(df_growth_EPX, ID %in% ID_alive_EPX)    
  
  df_growth_IMD <- subset(df_growth_tot, Molec %in%                      
                            c("Ctrl_1", "Ctrl_2", "IMD")) |>            
    mutate(Molecule = "IMD")                                              
  ID_alive_IMD <- subset(df_growth_IMD, t==28 & (!Status=="D"))$ID       
  df_growth_IMD_alive <- subset(df_growth_IMD, ID %in% ID_alive_IMD)      
  
  df_growth <- rbind(df_growth_EPX, df_growth_IMD)                      
  df_growth_alive <- rbind(df_growth_EPX_alive, df_growth_IMD_alive)       
  
  df_res <- list(
    df_expA_growth = df_growth,
    df_expA_growth_alive = df_growth_alive
  )
  
  return(df_res)
}

# Mixture - Jonker interaction models ----

# Fonction calculating the surface dose - response knowing the dose - response curves

CA_complete2 <- function(C_mat, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none", multicore = FALSE, mc.cores = 4) {
  param <- data.frame(Slopes, Ec50s)
  # try(cat(Slopes))
  # try(cat(Ec50s))

  # 1. CA model ----
  if (interact == "none") {
    a <- 0
    b <- 0
    f <- function(Y, x = x, Max = Max, param = param, a = a, b = b) {
      ecs <- param$Ec50s * ((Max - Y) / Y)^(1 / param$Slopes) # C1 corresponding to Y
      G <- (sum(x / ecs)) - 1
      return(abs(G))
    }
  }
  # 2. SA model ----
  else if (interact == "SA") {
    b <- 0
    f <- function(Y, x = x, Max = Max, param = param, a = a, b = b) {
      ecs <- param$Ec50s * ((Max - Y) / Y)^(1 / param$Slopes) # C1 corresponding to Y
      G <- (sum(x / ecs)) - exp(a * prod(x / param$Ec50s / (sum(x / param$Ec50s))))
      return(abs(G))
    }
  }
  # 3. DR model ----
  else if (interact == "DR") {
    if (length(b) != (dim(C_mat)[2] - 1)) {
      stop("length of b should be equal to the number of chemicals -1.")
    }
    f <- function(Y, x = x, Max = Max, param = param, a = a, b = b) {
      ecs <- param$Ec50s * ((Max - Y) / Y)^(1 / param$Slopes) # Cs corresponding to Y
      G <- (sum(x / ecs)) - exp((a + b %*% ((x / param$Ec50s / (sum(x / param$Ec50s))))[1:(dim(C_mat)[2] - 1)]) * prod(x / param$Ec50s / (sum(x / param$Ec50s))))
      return(abs(G))
    }
  }
  # 4. DR2 model ----
  else if (interact == "DR2") {
    if (length(b) != (dim(C_mat)[2] - 1)) {
      stop("length of b should be equal to the number of chemicals -1.")
    }
    f <- function(Y, x = x, Max = Max, param = param, a = a, b = b) {
      ecs <- param$Ec50s * ((Max - Y) / Y)^(1 / param$Slopes) # Cs corresponding to Y
      G <- (sum(x / ecs)) - exp((a + b %*% sin((x / param$Ec50s / (sum(x / param$Ec50s))) * 2 * pi)[1:(dim(C_mat)[2] - 1)]) * prod(x / param$Ec50s / (sum(x / param$Ec50s))))
      return(abs(G))
    }
  }
  # 5. DL model ----
  else if (interact == "DL") {
    f <- function(Y, x = x, Max = Max, param = param, a = a, b = b) {
      if (length(b) != (dim(C_mat)[2] - 1)) {
        stop("length of b should be equal to the number of chemicals -1.")
      }
      ecs <- param$Ec50s * ((Max - Y) / Y)^(1 / param$Slopes) # C1 corresponding to Y
      G <- (sum(x / ecs)) - exp(a * (1 - b * (sum(x / param$Ec50s))) * prod(x / param$Ec50s / (sum(x / param$Ec50s))))
      return(abs(G))
    }
  } else {
    stop("please specify interaction model")
  }

  Y_f <- function(x, Max = Max, param = param, a = a, b = b) {
    if (all(x == 0)) { # no chemical
      if (all(param$Slopes > 0)) {
        return(Max)
      }
      if (all(param$Slopes < 0)) {
        return(0)
      }
    }
    # Y<-optimize(f=function(Y) f(Y=Y, x=x, Max=Max, param=param, a=a, b=b), interval=c(0,Max*1.2))$minimum
    Y <- optimize(f = function(Y) f(Y = Y, x = x, Max = Max, param = param, a = a, b = b), interval = c(0, Max))$minimum

    return(Y)
  }
  if (multicore) {
    C_mat_list <- split(
      unique(C_mat), 
      cut(1:(dim(unique(C_mat))[1]), 
          breaks = dim(unique(C_mat))[1])
      )
    res_CA_unique <- unlist(
      mclapply(
        C_mat_list, 
        function(x) Y_f(x, Max = Max, param = param, a = a, b = b), 
        mc.cores = mc.cores
        )
      )
    res_CA <- rep(NA, dim(unique(C_mat))[1])
    for (i in 1:(dim(unique(C_mat))[1])) {
      res_CA[which((C_mat[, 1] == unique(C_mat)[i, 1]) & (C_mat[, 2] == unique(C_mat)[i, 2]))] <- res_CA_unique[i]
    }
  } else {
    res_CA_unique <- apply(
      unique(C_mat), 
      1, 
      function(x) Y_f(x, Max = Max, param = param, a = a, b = b)
      )
    res_CA <- rep(NA, dim(unique(C_mat))[1])
    for (i in 1:(dim(unique(C_mat))[1])) {
      res_CA[which((C_mat[, 1] == unique(C_mat)[i, 1]) & (C_mat[, 2] == unique(C_mat)[i, 2]))] <- res_CA_unique[i]
    }
  }
  if (any(is.na(res_CA))) {
    cat(param, a, b)
  }
  return(res_CA)
}

CA_complete2_fit_speed <- function(C_mat, Response, param = NULL, upper = NULL, lower = NULL, start = NULL, interact = "none", identical_slopes = FALSE, error_type = "Normal", iter = 500, multicore = FALSE, mc.cores = 4) {
  
  tol_drc <- 10^-8
  
  # 1. If dose-response curves known ----
  if (!is.null(param)) {
    if ((is.null(param$Max)) | (is.null(param$Slopes)) | (is.null(param$Ec50s))) {
      stop("param misspecification")
    }
    if (interact == "none") {
      stop("please specify interaction")
    }

    ## 1.1. SA model ----
    if (interact == "SA") {
      if (missing(upper)) {
        upper <- 20
      } # the intervals for a and b must be larger because they can compensate each other
      if (missing(lower)) {
        lower <- -20
      }

      # fit
      if (error_type == "Normal") {
        res_CA_optim <- optimize(
          f = function(x) {
            CA_complete2_RSS(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x,
              interact  = interact,
              multicore = multicore,
              mc.cores = mc.cores
            )
          },
          upper = upper,
          lower = lower
        )

        Res <- list(
          a = res_CA_optim$minimum,
          Error = res_CA_optim$objective
        )
        return(Res)
      } else if (error_type == "Poisson") {
        res_CA_optim <- optimize(
          f = function(x) {
            CA_complete2_Poisson(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x,
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          },
          upper = upper,
          lower = lower
        )

        Res <- list(
          a     = res_CA_optim$minimum,
          Error = res_CA_optim$objective
        )

        return(Res)
      } else if (error_type == "NB") {
        res_CA_optim <- optimize(
          f = function(x) {
            CA_complete2_NB(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x,
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          },
          upper = upper,
          lower = lower
        )
        
        Res <- list(
          a     = res_CA_optim$minimum,
          Error = res_CA_optim$objective
        )
        
        return(Res)
      } else {
        stop("Misspecification of the error model")
      }
    } # End SA model
  } # End dose-response curves known

  # 2. If dose-response curves not known ----
  if (is.null(param)) {
    require(dfoptim)
    
    print("estimation drc")

    if (any(C_mat == 0)) {
      mean_C_mat <- exp(apply(log(C_mat[-which(C_mat == 0, arr.ind = TRUE)[, 1], ]), 2, mean))
    } else {
      mean_C_mat <- exp(apply(log(C_mat), 2, mean))
    }

    ## 2.1. Different slopes ----
    if (!identical_slopes) {
      ### 2.1.1. CA model ----
      if (interact == "none") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 0, 0, mean_C_mat * 10)
        }

        if (missing(lower)) {
          lower <- c(0, -100, -100, 0, 0)
        }

        if (missing(start)) {
          start <- c(max(Response), -1, -1, mean_C_mat)
        }

        # fit
        if (error_type == "Normal") {
          print("nk common slopes Normal CA")
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )

          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          print("nk diff slopes Poisson CA")
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )

          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          # print("nmk diff slopes NB CA")
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            Error  = res_CA_nmk$value
          )
          return(Res)
        }else {
          stop("Misspecification of the error model")
        }
      } # End CA model

      ### 2.1.2. SA model ----
      if (interact == "SA") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 0, 0, mean_C_mat * 10, 20)
        } # the intervals for a and b must be larger beacuse they can compensate each other

        if (missing(lower)) {
          lower <- c(0, -100, -100, 0, 0, -20)
        }

        if (missing(start)) {
          start <- c(max(Response), -1, -1, mean_C_mat, 0)
        }

        # fit
        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                a         = x[6],
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )

          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            a      = res_CA_nmk$par[6],
            Error  = res_CA_nmk$value
          )

          return(Res)
        } else if (error_type == "Poisson") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                a         = x[6],
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )

          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            a      = res_CA_nmk$par[6],
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                a         = x[6],
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            a      = res_CA_nmk$par[6],
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End SA model
    } # End common slopes

    ## 2.2. Common slopes ----
    else {
      ### 2.2.1. CA model ----
      if (interact == "none") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10)
        }

        if (missing(lower)) {
          lower <- c(0, 0, 0, 0)
        }

        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat)
        }

        # fit
        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[2]),
                Ec50s     = c(x[3], x[4]),
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )

          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]),
            Ec50s  = c(res_CA_nmk$par[3:4]),
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          # print("nk common slopes Poisson CA")
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[2]),
                Ec50s     = c(x[3], x[4]),
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )

          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]),
            Ec50s  = c(res_CA_nmk$par[3:4]),
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          # print("nk common slopes NB CA")
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[2]),
                Ec50s     = c(x[3], x[4]),
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )
          
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]),
            Ec50s  = c(res_CA_nmk$par[3:4]),
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End CA model

      ### 2.2.2. SA model ----
      if (interact == "SA") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10, 20)
        } # the intervals for a and b must be larger beacuse they can compensate each other

        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, -20)
        }

        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat, 0)
        }

        # fit

        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[2]),
                Ec50s     = c(x[3], x[4]),
                a         = x[5],
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )

          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]),
            Ec50s  = c(res_CA_nmk$par[3:4]),
            a      = res_CA_nmk$par[5],
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[2]),
                Ec50s     = c(x[3], x[4]),
                a         = x[5],
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )

          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]),
            Ec50s  = c(res_CA_nmk$par[3:4]),
            a      = res_CA_nmk$par[5],
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[2]),
                Ec50s     = c(x[3], x[4]),
                a         = x[5],
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )
          
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]),
            Ec50s  = c(res_CA_nmk$par[3:4]),
            a      = res_CA_nmk$par[5],
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End SA model
    } # End different slopes
  } # End dose-response curves not known
} # End.

CA_complete2_RSS <- function(C_mat, Response, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none", multicore = FALSE, mc.cores = 4) {
  
  if (length(Response) != (dim(C_mat)[1])) {
    stop("Should be as many responses as there are conditions")
  }

  res_CA <- CA_complete2(
    C_mat     = C_mat,
    Max       = Max,
    Slopes    = Slopes,
    Ec50s     = Ec50s,
    a         = a,
    b         = b,
    interact  = interact,
    multicore = multicore,
    mc.cores  = mc.cores
  )
  res_CA_RSS <- sum((res_CA - Response)^2)

  return(res_CA_RSS)
}


CA_complete2_Poisson <- function(C_mat, Response, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none", multicore = FALSE, mc.cores = 4) {
  
  if (length(Response) != (dim(C_mat)[1])) {
    stop("Should be as many responses as there are conditions")
  }

  res_CA <- CA_complete2(
    C_mat     = C_mat,
    Max       = Max,
    Slopes    = Slopes,
    Ec50s     = Ec50s,
    a         = a,
    b         = b,
    interact  = interact,
    multicore = multicore,
    mc.cores  = mc.cores
  )

  # Pour stabilité numérique, éviter log(0)
  if (any(res_CA <= 0)) {
    return(-Inf)
  }
  # res_CA_Poisson <- -sum(Response * log(res_CA) - res_CA - lfactorial(Response)) # - Loglikelihood
  res_CA_Poisson <- -sum(dpois(Response, lambda = res_CA, log = TRUE))
  
  return(res_CA_Poisson)
}

CA_complete2_QuasiPoisson <- function(C_mat, Response, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none", multicore = FALSE, mc.cores = 4) {
  
  if (length(Response) != (dim(C_mat)[1])) {
    stop("Should be as many responses as there are conditions")
  }
  
  res_CA <- CA_complete2(
    C_mat     = C_mat,
    Max       = Max,
    Slopes    = Slopes,
    Ec50s     = Ec50s,
    a         = a,
    b         = b,
    interact  = interact,
    multicore = multicore,
    mc.cores  = mc.cores
  )
  
  # Pour stabilité numérique, éviter log(0)
  if (any(res_CA <= 0)) {
    return(-Inf)
  }
  
  # Number of parameters
  if (interact == "none"){
    p <- 0
  } else if (interact == "SA") {
    p <- 1
  } else if (interact %in% c("DR", "DL")) {
    p <- 2
  }
  
  term1 <- ifelse(Response == 0, 0, Response * log(Response / res_CA))
  term2 <- Response - res_CA
  deviance_Poisson <- 2 * sum(term1 - term2)
  
  df_residual <- length(Response) - p
  dev <- deviance_quasipoisson(Response, res_CA)
  phi <- dev / df_residual
  
  # Pseudo Log-Vraissemblance
  res_CA_QuasiPoisson <- -sum(Response * log(res_CA) - res_CA - lfactorial(Response)) / phi # - Loglikelihood
  
  return(res_CA_QuasiPoisson)
}

CA_complete2_NB <- function(C_mat, Response, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none", multicore = FALSE, mc.cores = 4) {
  
  if (length(Response) != (dim(C_mat)[1])) {
    stop("Should be as many responses as there are conditions")
  }
  
  res_CA <- CA_complete2(
    C_mat     = C_mat,
    Max       = Max,
    Slopes    = Slopes,
    Ec50s     = Ec50s,
    a         = a,
    b         = b,
    interact  = interact,
    multicore = multicore,
    mc.cores  = mc.cores
  )
  
  # Pour stabilité numérique, éviter log(0)
  if (any(res_CA <= 0)) {
    return(-Inf)
  }
  
  neg_loglik_for_theta <- function(log_theta, y, mu) {
    theta <- exp(log_theta)  # contrainte >0
    -sum(dnbinom(y, size = theta, mu = res_CA, log = TRUE))
  }
  res <- optim(par = log(1), fn = neg_loglik_for_theta, y = Response, mu = res_CA,
               method = "Brent", lower = log(1e-8), upper = log(1e6)) 
  
  est_theta <- exp(res$par)
  print(est_theta)
  
  res_CA_NB <- -sum(dnbinom(Response, size = est_theta, mu = res_CA, log = TRUE))[1]
  
  # Log-Vraissemblance
  
  return(res_CA_NB)
}

CA_complete_fit_speed <- function(C_mat, Response, param = NULL, upper = NULL, lower = NULL, start = NULL, interact = "none", identical_slopes = FALSE, error_type = "Normal", iter = 500, multicore = FALSE, mc.cores = 4) {
  
  tol_drc <- 10^-8
  
  # 1. Dose-response curves known ----
  if (!is.null(param)) {
    if ((is.null(param$Max)) | (is.null(param$Slopes)) | (is.null(param$Ec50s))) {
      stop("param misspecification")
    }
    if (interact == "none") {
      stop("please specify interaction")
    }

    ## 1.1. SA model ----
    if (interact == "SA") {
      if (missing(upper)) {
        upper <- 20
      } # the intervals for a and b must be larger because they can compensate each other
      if (missing(lower)) {
        lower <- -20
      }

      if (error_type == "Normal") {
        res_CA_optim <- optimize(
          f = function(x) {
            CA_complete2_RSS(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x,
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          },
          upper = upper,
          lower = lower
        )
        Res <- list(
          a     = res_CA_optim$minimum,
          Error = res_CA_optim$objective
        )
        return(Res)
      } else if (error_type == "Poisson") {
        res_CA_optim <- optimize(
          f = function(x) {
            CA_complete2_Poisson(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x,
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          },
          upper = upper,
          lower = lower
        )
        Res <- list(
          a     = res_CA_optim$minimum,
          Error = res_CA_optim$objective
        )
        return(Res)
      } else if (error_type == "NB") {
        res_CA_optim <- optimize(
          f = function(x) {
            CA_complete2_NB(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x,
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          },
          upper = upper,
          lower = lower
        )
        Res <- list(
          a     = res_CA_optim$minimum,
          Error = res_CA_optim$objective
        )
        return(Res)
      } else {
        stop("Misspecification of the error model")
      }
    } # End SA

    ## 1.2. DR model ----
    if (interact == "DR") {
      if (error_type == "Normal") {
        res_CA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            CA_complete2_RSS(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          }
        )
        cat(res_CA_optim$convergence)
        Res <- list(
          a     = res_CA_optim$par[1],
          b     = res_CA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_CA_optim$value
        )
        return(Res)
      } else if (error_type == "Poisson") {
        res_CA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            CA_complete2_Poisson(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          }
        )
        cat(res_CA_optim$convergence)

        Res <- list(
          a     = res_CA_optim$par[1],
          b     = res_CA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_CA_optim$value
        )
        return(Res)
      } else if (error_type == "NB") {
        res_CA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            CA_complete2_NB(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          }
        )
        cat(res_CA_optim$convergence)
        
        Res <- list(
          a     = res_CA_optim$par[1],
          b     = res_CA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_CA_optim$value
        )
        return(Res)
      } else {
        stop("Misspecification of the error model")
      }
    } # End DR

    ## 1.3. DR2 model ----
    if (interact == "DR2") {
      if (error_type == "Normal") {
        res_CA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            CA_complete2_RSS(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          }
        )
        cat(res_CA_optim$convergence)

        Res <- list(
          a     = res_CA_optim$par[1],
          b     = res_CA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_CA_optim$value
        )
        return(Res)
      } else if (error_type == "Poisson") {
        res_CA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            CA_complete2_Poisson(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          }
        )
        cat(res_CA_optim$convergence)
        Res <- list(
          a     = res_CA_optim$par[1],
          b     = res_CA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_CA_optim$value
        )
        return(Res)
      } else if (error_type == "NB") {
        res_CA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            CA_complete2_NB(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          }
        )
        cat(res_CA_optim$convergence)
        Res <- list(
          a     = res_CA_optim$par[1],
          b     = res_CA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_CA_optim$value
        )
        return(Res)
      } else {
        stop("Misspecification of the error model")
      }
    } # End DR2

    ## 1.4. DL model ----
    if (interact == "DL") {
      if (error_type == "Normal") {
        res_CA_optim <- optim(
          par = c(0, 0),
          fn = function(x) {
            CA_complete2_RSS(
              C_mat = C_mat,
              Response = Response,
              Max = param$Max,
              Slopes = param$Slopes,
              Ec50s = param$Ec50s,
              a = x[1], b = x[2],
              interact = interact,
              multicore = multicore,
              mc.cores = mc.cores
            )
          }
        )
        cat(res_CA_optim$convergence)
        Res <- list(
          a     = res_CA_optim$par[1],
          b     = res_CA_optim$par[2],
          Error = res_CA_optim$value
        )
        return(Res)
      } else if (error_type == "Poisson") {
        res_CA_optim <- optim(
          par = c(0, 0),
          fn = function(x) {
            CA_complete2_Poisson(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2],
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          }
        )
        cat(res_CA_optim$convergence)
        Res <- list(
          a     = res_CA_optim$par[1],
          b     = res_CA_optim$par[2],
          Error = res_CA_optim$value
        )
        return(Res)
      } else if (error_type == "NB") {
        res_CA_optim <- optim(
          par = c(0, 0),
          fn = function(x) {
            CA_complete2_NB(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2],
              interact  = interact,
              multicore = multicore,
              mc.cores  = mc.cores
            )
          }
        )
        cat(res_CA_optim$convergence)
        Res <- list(
          a     = res_CA_optim$par[1],
          b     = res_CA_optim$par[2],
          Error = res_CA_optim$value
        )
        return(Res)
      } else {
        stop("Misspecification of the error model")
      }
    } # End DL
  } # End dose-response curves known

  # 2. Dose-response curves not known ----
  if (is.null(param)) {
    # require(DEoptim)
    require(dfoptim)
    if (any(C_mat == 0)) {
      mean_C_mat <- exp(apply(log(C_mat[-which(C_mat == 0, arr.ind = TRUE)[, 1], ]), 2, mean))
    } else {
      mean_C_mat <- exp(apply(log(C_mat), 2, mean))
    }

    ## 2.1. Different slopes ----
    if (!identical_slopes) {
      ### 2.1.1. CA model ----
      if (interact == "none") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, 100, mean_C_mat * 10)
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, 0)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, 1, mean_C_mat)
        }

        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      }

      ### 2.1.2. SA model ----
      if (interact == "SA") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, 100, mean_C_mat * 10, 20) # the intervals for a and b must be larger beacuse they can compensate each other
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, 0, -20)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, 1, mean_C_mat, 0)
        }

        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                a         = x[6],
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]),
            a      = res_CA_nmk$par[6],
            Error  = res_CA_nmk$value
          )
          return(Res)
          
        } else if (error_type == "Poisson") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]), 
            Ec50s  = c(res_CA_nmk$par[4:5]), 
            a      = res_CA_nmk$par[6], 
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_CA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2:3]), 
            Ec50s  = c(res_CA_nmk$par[4:5]), 
            a      = res_CA_nmk$par[6], 
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End SA 
      
      ### 2.1.3. DR model ----
      if (interact == "DR") {
        
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, 100, mean_C_mat * 10, 20, rep(30, (dim(C_mat)[2] - 1)))
        } # the intervals for a and b must be larger beacuse they can compensate each other
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, 0, -20, rep(-30, (dim(C_mat)[2] - 1)))
        }
        if (missing(start)) {
          start <- c(max(Response), 1, 1, mean_C_mat, 0, rep(0, (dim(C_mat)[2] - 1)))
        }

        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], 
                b         = x[7:length(x)], 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
                )
              }, 
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2:3]),
            Ec50s  = c(res_CA_nmk$par[4:5]), 
            a      = res_CA_nmk$par[6], 
            b      = res_CA_nmk$par[7:(7 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "Poisson") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]), 
                a         = x[6],
                b         = x[7:length(x)], 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
                )
              },
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2:3]), 
            Ec50s  = c(res_CA_nmk$par[4:5]), 
            a      = res_CA_nmk$par[6], 
            b      = res_CA_nmk$par[7:(7 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "NB") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]), 
                a         = x[6],
                b         = x[7:length(x)], 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
              )
            },
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2:3]), 
            Ec50s  = c(res_CA_nmk$par[4:5]), 
            a      = res_CA_nmk$par[6], 
            b      = res_CA_nmk$par[7:(7 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End DR 
      
      ### 2.1.4. DL model ----
      if (interact == "DL") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, 100, mean_C_mat * 10, 20, 30)
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, 0, -20, -30)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, 1, mean_C_mat, 0, 0)
        }

        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], b = x[7], 
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
                )
              }, 
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2:3]), 
            Ec50s  = c(res_CA_nmk$par[4:5]), 
            a      = res_CA_nmk$par[6], 
            b      = res_CA_nmk$par[7], 
            Error  = res_CA_nmk$value
            )
          return(Res)
          
        } else if (error_type == "Poisson") {
          
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], b = x[7],
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
                )
              }, 
              upper   = upper, 
              lower   = lower, 
              control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2:3]), 
            Ec50s  = c(res_CA_nmk$par[4:5]), 
            a      = res_CA_nmk$par[6], 
            b      = res_CA_nmk$par[7], 
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "NB") {
          
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], b = x[7],
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
              )
            }, 
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2:3]), 
            Ec50s  = c(res_CA_nmk$par[4:5]), 
            a      = res_CA_nmk$par[6], 
            b      = res_CA_nmk$par[7], 
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End DL
    } 
    ## 2.2. Common slopes ----
    else {
      
      ### 2.2.1. CA model ----
      if (interact == "none") {
        
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10)
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat)
        }

        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3], x[4]), 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
                )
              }, 
            upper = upper, 
            lower = lower,
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "Poisson") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3], x[4]), 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
                )
              }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "NB") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3], x[4]), 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End CA
      
      ### 2.2.2. SA model ----
      if (interact == "SA") {
        
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10, 20)
        } # the intervals for a and b must be larger beacuse they can compensate each other
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, -20)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat, 0)
        }

        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], Slopes = c(x[2], x[2]), 
                Ec50s     = c(x[3], x[4]), a = x[5], 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
                )
              }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]),
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            a      = res_CA_nmk$par[5], 
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "Poisson") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat = C_mat, 
                Response = Response, 
                Max = x[1], 
                Slopes = c(x[2], x[2]), 
                Ec50s = c(x[3], x[4]), 
                a = x[5], 
                interact = interact,
                multicore = multicore,
                mc.cores = mc.cores
                )
              }, 
            upper = upper, 
            lower = lower,
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s = c(res_CA_nmk$par[3:4]), 
            a = res_CA_nmk$par[5],
            Error = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "NB") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_NB(
                C_mat = C_mat, 
                Response = Response, 
                Max = x[1], 
                Slopes = c(x[2], x[2]), 
                Ec50s = c(x[3], x[4]), 
                a = x[5], 
                interact = interact,
                multicore = multicore,
                mc.cores = mc.cores
              )
            }, 
            upper = upper, 
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s = c(res_CA_nmk$par[3:4]), 
            a = res_CA_nmk$par[5],
            Error = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End SA
      
      ### 2.2.3. DR model ----
      if (interact == "DR") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10, 20, rep(30, (dim(C_mat)[2] - 1)))
        } # the intervals for a and b must be larger beacuse they can compensate each other
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, -20, rep(-30, (dim(C_mat)[2] - 1)))
        }
        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat, 0, 0)
        }

        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[2]),
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6:length(x)],
                interact  = interact,
                multicore = multicore, 
                mc.cores  = mc.cores
                )
              }, 
            upper = upper,
            lower = lower, 
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1],
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            a      = res_CA_nmk$par[5],
            b      = res_CA_nmk$par[6:(6 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "Poisson") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6:length(x)], 
                interact  = interact, 
                multicore = multicore,
                mc.cores  = mc.cores
                )
              }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            a      = res_CA_nmk$par[5],
            b      = res_CA_nmk$par[6:(6 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "NB") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6:length(x)], 
                interact  = interact, 
                multicore = multicore,
                mc.cores  = mc.cores
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            a      = res_CA_nmk$par[5],
            b      = res_CA_nmk$par[6:(6 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      }
      
      ### 2.2.4. DL model ----
      if (interact == "DL") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10, 20, 30)
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, -20, -30)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat, 0, 0)
        }

        if (error_type == "Normal") {
          res_CA_nmk <- nmkb(
            par = start, fn = function(x) {
              CA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]),
                a         = x[5], 
                b         = x[6], 
                interact  = interact,
                multicore = multicore,
                mc.cores  = mc.cores
                )
              }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]),
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            a      = res_CA_nmk$par[5], 
            b      = res_CA_nmk$par[6],
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "Poisson") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6], 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
                )
              }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
            )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            a      = res_CA_nmk$par[5],
            b      = res_CA_nmk$par[6], 
            Error  = res_CA_nmk$value
            )
          return(Res)
        } else if (error_type == "NB") {
          res_CA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              CA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6], 
                interact  = interact, 
                multicore = multicore, 
                mc.cores  = mc.cores
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_CA_nmk$par[1], 
            Slopes = c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
            Ec50s  = c(res_CA_nmk$par[3:4]), 
            a      = res_CA_nmk$par[5],
            b      = res_CA_nmk$par[6], 
            Error  = res_CA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End DL
    } # End common slopes
  } # End Dose-response curves not known
} # End.






IA_complete2 <- function(C_mat, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none") {
  param <- data.frame(Slopes, Ec50s)
  # try(cat(Slopes))
  # try(cat(Ec50s))
  
  # 1. IA model ----
  
  if (interact == "none"){
    
    Y <- c()
    
    for (i in 1:nrow(C_mat)){

      yA_i <- 1/(1+exp(as.numeric(Slopes[1])*(log(C_mat[i,1])-log(as.numeric(Ec50s[1])))))
      yB_i <- 1/(1+exp(as.numeric(Slopes[2])*(log(C_mat[i,2])-log(as.numeric(Ec50s[2])))))
      
      Y_i <- as.numeric(Max[1]) * (yA_i*yB_i)
      Y <- c(Y, Y_i)
    }
  
  # 2. SA model ----
  
  } else if (interact == "SA") {
    
    Y <- c()
    
    for (i in 1:nrow(C_mat)){
      yA_i <- 1/(1+exp(as.numeric(Slopes[1])*(log(C_mat[i,1])-log(as.numeric(Ec50s[1])))))
      yB_i <- 1/(1+exp(as.numeric(Slopes[2])*(log(C_mat[i,2])-log(as.numeric(Ec50s[2])))))
      
      TUA_i <- C_mat[i,1]/as.numeric(Ec50s[1])
      TUB_i <- C_mat[i,2]/as.numeric(Ec50s[2])
      
      zA_i <- TUA_i/(TUA_i + TUB_i)
      zB_i <- TUB_i/(TUA_i + TUB_i)
      
      G <- a * (zA_i*zB_i)
      
      Y_i <- as.numeric(Max[1]) * pnorm(qnorm(yA_i*yB_i) + G)
      Y <- c(Y, Y_i)
    }
    
  # 3. DR model ----
  
  } else if (interact == "DR"){
    
    Y <- c()
    
    for (i in 1:nrow(C_mat)){
      yA_i <- 1/(1+exp(as.numeric(Slopes[1])*(log(C_mat[i,1])-log(as.numeric(Ec50s[1])))))
      yB_i <- 1/(1+exp(as.numeric(Slopes[2])*(log(C_mat[i,2])-log(as.numeric(Ec50s[2])))))
      
      TUA_i <- C_mat[i,1]/as.numeric(Ec50s[1])
      TUB_i <- C_mat[i,2]/as.numeric(Ec50s[2])
      
      zA_i <- TUA_i/(TUA_i + TUB_i)
      zB_i <- TUB_i/(TUA_i + TUB_i)
      
      G <- (a+b*zA_i) * (zA_i*zB_i)
      
      Y_i <- as.numeric(Max[1]) * pnorm(qnorm(yA_i*yB_i) + G)
      Y <- c(Y, Y_i)
    }
  
  
  # 4. DL model ----
  
  } else if (interact == "DL"){
    
    Y <- c()
    
    for (i in 1:nrow(C_mat)){
      yA_i <- 1/(1+exp(as.numeric(Slopes[1])*(log(C_mat[i,1])-log(as.numeric(Ec50s[1])))))
      yB_i <- 1/(1+exp(as.numeric(Slopes[2])*(log(C_mat[i,2])-log(as.numeric(Ec50s[2])))))
      
      TUA_i <- C_mat[i,1]/as.numeric(Ec50s[1])
      TUB_i <- C_mat[i,2]/as.numeric(Ec50s[2])
      
      zA_i <- TUA_i/(TUA_i + TUB_i)
      zB_i <- TUB_i/(TUA_i + TUB_i)
      
      G <- a*(1-b*(yA_i*yB_i)) * zA_i*zB_i
      
      Y_i <- as.numeric(Max[1]) * pnorm(qnorm(yA_i*yB_i) + G)
      Y <- c(Y, Y_i)
    }
  }
  
  res_IA <- Y

  return(res_IA)
}


IA_complete2_RSS <- function(C_mat, Response, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none") {
  
  if (length(Response) != (dim(C_mat)[1])) {
    stop("Should be as many responses as there are conditions")
  }
  
  res_IA <- IA_complete2(
    C_mat     = C_mat,
    Max       = Max,
    Slopes    = Slopes,
    Ec50s     = Ec50s,
    a         = a,
    b         = b,
    interact  = interact,
    multicore = multicore,
    mc.cores  = mc.cores
  )
  res_IA_RSS <- sum((res_IA - Response)^2)
  
  return(res_IA_RSS)
}


IA_complete2_Poisson <- function(C_mat, Response, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none") {
  
  if (length(Response) != (dim(C_mat)[1])) {
    stop("Should be as many responses as there are conditions")
  }
  
  res_IA <- IA_complete2(
    C_mat     = C_mat,
    Max       = Max,
    Slopes    = Slopes,
    Ec50s     = Ec50s,
    a         = a,
    b         = b,
    interact  = interact
  )
  
  # Pour stabilité numérique, éviter log(0)
  if (any(res_IA <= 0)) {
    return(-Inf)
  }
  # res_IA_Poisson <- -sum(Response * log(res_IA) - res_IA - lfactorial(Response)) # - Loglikelihood
  
  res_IA_Poisson <- -sum(dpois(Response, lambda = res_IA, log = TRUE)) # - Loglikelihood
  
  return(res_IA_Poisson)
}

IA_complete2_QuasiPoisson <- function(C_mat, Response, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none") {
  
  if (length(Response) != (dim(C_mat)[1])) {
    stop("Should be as many responses as there are conditions")
  }
  
  res_IA <- IA_complete2(
    C_mat     = C_mat,
    Max       = Max,
    Slopes    = Slopes,
    Ec50s     = Ec50s,
    a         = a,
    b         = b,
    interact  = interact
  )
  
  # Pour stabilité numérique, éviter log(0)
  if (any(res_IA <= 0)) {
    return(-Inf)
  }
  
  # Number of parameters
  if (interact == "none"){
    p <- 0
  } else if (interact == "SA") {
    p <- 1
  } else if (interact %in% c("DR", "DL")) {
    p <- 2
  }
  
  term1 <- ifelse(Response == 0, 0, Response * log(Response / res_IA))
  term2 <- Response - res_IA
  deviance_Poisson <- 2 * sum(term1 - term2)
  
  df_residual <- length(Response) - p
  dev <- deviance_quasipoisson(Response, res_IA)
  phi <- dev / df_residual
  
  # Pseudo Log-Vraissemblance
  res_IA_QuasiPoisson <- -sum(Response * log(res_IA) - res_IA - lfactorial(Response)) / phi # - Loglikelihood
  
  return(res_IA_QuasiPoisson)
}

IA_complete2_NB <- function(C_mat, Response, Max, Slopes, Ec50s, a = 0, b = 0, interact = "none") {
  
  if (length(Response) != (dim(C_mat)[1])) {
    stop("Should be as many responses as there are conditions")
  }
  
  res_IA <- IA_complete2(
    C_mat     = C_mat,
    Max       = Max,
    Slopes    = Slopes,
    Ec50s     = Ec50s,
    a         = a,
    b         = b,
    interact  = interact
  )
  
  # Pour stabilité numérique, éviter log(0)
  if (any(res_IA <= 0)) {
    return(-Inf)
  }
  
  neg_loglik_for_theta <- function(log_theta, y, mu) {
    theta <- exp(log_theta)  # contrainte >0
    -sum(dnbinom(y, size = theta, mu = res_IA, log = TRUE))
  }
  res <- optim(par = log(1), fn = neg_loglik_for_theta, y = Response, mu = res_IA,
               method = "Brent", lower = log(1e-8), upper = log(1e6)) 
  
  est_theta <- exp(res$par)
  print(est_theta)
  # Log-Vraissemblance
  res_IA_NB <- -sum(dnbinom(Response, size = est_theta, mu = res_IA, log = TRUE)) # - Loglikelihood
  
  return(res_IA_NB)
}

IA_complete_fit_speed <- function(C_mat, Response, param = NULL, upper = NULL, lower = NULL, start = NULL, interact = "none", identical_slopes = FALSE, error_type = "Normal", iter = 500, multicore = FALSE, mc.cores = 4) {
  
  tol_drc <- 10^-8
  
  # 1. Dose-response curves known ----
  if (!is.null(param)) {
    if ((is.null(param$Max)) | (is.null(param$Slopes)) | (is.null(param$Ec50s))) {
      stop("param misspecification")
    }
    if (interact == "none") {
      stop("please specify interaction")
    }
    
    ## 1.1. SA model ----
    if (interact == "SA") {
      if (missing(upper)) {
        upper <- 20
      } # the intervals for a and b must be larger because they can compensate each other
      if (missing(lower)) {
        lower <- -20
      }
      
      if (error_type == "Normal") {
        res_IA_optim <- optimize(
          f = function(x) {
            IA_complete2_RSS(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x,
              interact  = interact
            )
          },
          upper = upper,
          lower = lower
        )
        Res <- list(
          a     = res_IA_optim$minimum,
          Error = res_IA_optim$objective
        )
        return(Res)
      } else if (error_type == "Poisson") {
        res_IA_optim <- optimize(
          f = function(x) {
            IA_complete2_Poisson(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x,
              interact  = interact
            )
          },
          upper = upper,
          lower = lower
        )
        Res <- list(
          a     = res_IA_optim$minimum,
          Error = res_IA_optim$objective
        )
        return(Res)
      } else if (error_type == "NB") {
        res_IA_optim <- optimize(
          f = function(x) {
            IA_complete2_NB(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x,
              interact  = interact
            )
          },
          upper = upper,
          lower = lower
        )
        Res <- list(
          a     = res_IA_optim$minimum,
          Error = res_IA_optim$objective
        )
        return(Res)
      } else {
        stop("Misspecification of the error model")
      }
    } # End SA
    
    ## 1.2. DR model ----
    if (interact == "DR") {
      if (error_type == "Normal") {
        res_IA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            IA_complete2_RSS(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact
            )
          }
        )
        cat(res_IA_optim$convergence)
        Res <- list(
          a     = res_IA_optim$par[1],
          b     = res_IA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_IA_optim$value
        )
        return(Res)
      } else if (error_type == "Poisson") {
        res_IA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            IA_complete2_Poisson(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact
            )
          }
        )
        cat(res_IA_optim$convergence)
        
        Res <- list(
          a     = res_IA_optim$par[1],
          b     = res_IA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_IA_optim$value
        )
        return(Res)
      } else if (error_type == "NB") {
        res_IA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            IA_complete2_NB(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact
            )
          }
        )
        cat(res_IA_optim$convergence)
        
        Res <- list(
          a     = res_IA_optim$par[1],
          b     = res_IA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_IA_optim$value
        )
        return(Res)
      } else {
        stop("Misspecification of the error model")
      }
    } # End DR
    
    ## 1.3. DR2 model ----
    if (interact == "DR2") {
      if (error_type == "Normal") {
        res_IA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            IA_complete2_RSS(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact
            )
          }
        )
        cat(res_IA_optim$convergence)
        
        Res <- list(
          a     = res_IA_optim$par[1],
          b     = res_IA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_IA_optim$value
        )
        return(Res)
      } else if (error_type == "Poisson") {
        res_IA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            IA_complete2_Poisson(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact
            )
          }
        )
        cat(res_IA_optim$convergence)
        Res <- list(
          a     = res_IA_optim$par[1],
          b     = res_IA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_IA_optim$value
        )
        return(Res)
      } else if (error_type == "NB") {
        res_IA_optim <- optim(
          par = rep(0, dim(C_mat)[2]),
          fn = function(x) {
            IA_complete2_NB(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2:length(x)],
              interact  = interact
            )
          }
        )
        cat(res_IA_optim$convergence)
        Res <- list(
          a     = res_IA_optim$par[1],
          b     = res_IA_optim$par[2:(2 + (dim(C_mat)[2] - 1) - 1)],
          Error = res_IA_optim$value
        )
        return(Res)
      } else {
        stop("Misspecification of the error model")
      }
    } # End DR2
    
    ## 1.4. DL model ----
    if (interact == "DL") {
      if (error_type == "Normal") {
        res_IA_optim <- optim(
          par = c(0, 0),
          fn = function(x) {
            IA_complete2_RSS(
              C_mat = C_mat,
              Response = Response,
              Max = param$Max,
              Slopes = param$Slopes,
              Ec50s = param$Ec50s,
              a = x[1], b = x[2],
              interact = interact
            )
          }
        )
        cat(res_IA_optim$convergence)
        Res <- list(
          a     = res_IA_optim$par[1],
          b     = res_IA_optim$par[2],
          Error = res_IA_optim$value
        )
        return(Res)
      } else if (error_type == "Poisson") {
        res_IA_optim <- optim(
          par = c(0, 0),
          fn = function(x) {
            IA_complete2_Poisson(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2],
              interact  = interact
            )
          }
        )
        cat(res_IA_optim$convergence)
        Res <- list(
          a     = res_IA_optim$par[1],
          b     = res_IA_optim$par[2],
          Error = res_IA_optim$value
        )
        return(Res)
      } else if (error_type == "NB") {
        res_IA_optim <- optim(
          par = c(0, 0),
          fn = function(x) {
            IA_complete2_NB(
              C_mat     = C_mat,
              Response  = Response,
              Max       = param$Max,
              Slopes    = param$Slopes,
              Ec50s     = param$Ec50s,
              a         = x[1],
              b         = x[2],
              interact  = interact
            )
          }
        )
        cat(res_IA_optim$convergence)
        Res <- list(
          a     = res_IA_optim$par[1],
          b     = res_IA_optim$par[2],
          Error = res_IA_optim$value
        )
        return(Res)
      } else {
        stop("Misspecification of the error model")
      }
    } # End DL
  } # End dose-response curves known
  
  # 2. Dose-response curves not known ----
  if (is.null(param)) {
    # require(DEoptim)
    require(dfoptim)
    if (any(C_mat == 0)) {
      mean_C_mat <- exp(apply(log(C_mat[-which(C_mat == 0, arr.ind = TRUE)[, 1], ]), 2, mean))
    } else {
      mean_C_mat <- exp(apply(log(C_mat), 2, mean))
    }
    
    ## 2.1. Different slopes ----
    if (!identical_slopes) {
      ### 2.1.1. IA model ----
      if (interact == "none") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, 100, mean_C_mat * 10)
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, 0)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, 1, mean_C_mat)
        }
        
        if (error_type == "Normal") {
          res_IA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              IA_complete2_RSS(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                interact  = interact
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1],
            Slopes = c(res_IA_nmk$par[2:3]),
            Ec50s  = c(res_IA_nmk$par[4:5]),
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          res_IA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              IA_complete2_Poisson(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                interact  = interact
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1],
            Slopes = c(res_IA_nmk$par[2:3]),
            Ec50s  = c(res_IA_nmk$par[4:5]),
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_IA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              IA_complete2_NB(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                interact  = interact
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1],
            Slopes = c(res_IA_nmk$par[2:3]),
            Ec50s  = c(res_IA_nmk$par[4:5]),
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      }
      
      ### 2.1.2. SA model ----
      if (interact == "SA") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, 100, mean_C_mat * 10, 20) # the intervals for a and b must be larger beacuse they can compensate each other
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, 0, -20)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, 1, mean_C_mat, 0)
        }
        
        if (error_type == "Normal") {
          res_IA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              IA_complete2_RSS(
                C_mat     = C_mat,
                Response  = Response,
                Max       = x[1],
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]),
                a         = x[6],
                interact  = interact
              )
            },
            upper = upper,
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1],
            Slopes = c(res_IA_nmk$par[2:3]),
            Ec50s  = c(res_IA_nmk$par[4:5]),
            a      = res_IA_nmk$par[6],
            Error  = res_IA_nmk$value
          )
          return(Res)
          
        } else if (error_type == "Poisson") {
          res_IA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              IA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], 
                interact  = interact
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1],
            Slopes = c(res_IA_nmk$par[2:3]), 
            Ec50s  = c(res_IA_nmk$par[4:5]), 
            a      = res_IA_nmk$par[6], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_IA_nmk <- nmkb(
            par = start,
            fn = function(x) {
              IA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], 
                interact  = interact
              )
            },
            upper   = upper,
            lower   = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1],
            Slopes = c(res_IA_nmk$par[2:3]), 
            Ec50s  = c(res_IA_nmk$par[4:5]), 
            a      = res_IA_nmk$par[6], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End SA 
      
      ### 2.1.3. DR model ----
      if (interact == "DR") {
        
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, 100, mean_C_mat * 10, 20, rep(30, (dim(C_mat)[2] - 1)))
        } # the intervals for a and b must be larger beacuse they can compensate each other
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, 0, -20, rep(-30, (dim(C_mat)[2] - 1)))
        }
        if (missing(start)) {
          start <- c(max(Response), 1, 1, mean_C_mat, 0, rep(0, (dim(C_mat)[2] - 1)))
        }
        
        if (error_type == "Normal") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], 
                b         = x[7:length(x)], 
                interact  = interact
              )
            }, 
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2:3]),
            Ec50s  = c(res_IA_nmk$par[4:5]), 
            a      = res_IA_nmk$par[6], 
            b      = res_IA_nmk$par[7:(7 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]), 
                a         = x[6],
                b         = x[7:length(x)], 
                interact  = interact
              )
            },
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2:3]), 
            Ec50s  = c(res_IA_nmk$par[4:5]), 
            a      = res_IA_nmk$par[6], 
            b      = res_IA_nmk$par[7:(7 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]),
                Ec50s     = c(x[4], x[5]), 
                a         = x[6],
                b         = x[7:length(x)], 
                interact  = interact
              )
            },
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2:3]), 
            Ec50s  = c(res_IA_nmk$par[4:5]), 
            a      = res_IA_nmk$par[6], 
            b      = res_IA_nmk$par[7:(7 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End DR 
      
      ### 2.1.4. DL model ----
      if (interact == "DL") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, 100, mean_C_mat * 10, 20, 30)
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, 0, -20, -30)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, 1, mean_C_mat, 0, 0)
        }
        
        if (error_type == "Normal") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], b = x[7], 
                interact  = interact
              )
            }, 
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2:3]), 
            Ec50s  = c(res_IA_nmk$par[4:5]), 
            a      = res_IA_nmk$par[6], 
            b      = res_IA_nmk$par[7], 
            Error  = res_IA_nmk$value
          )
          return(Res)
          
        } else if (error_type == "Poisson") {
          
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], b = x[7],
                interact  = interact
              )
            }, 
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2:3]), 
            Ec50s  = c(res_IA_nmk$par[4:5]), 
            a      = res_IA_nmk$par[6], 
            b      = res_IA_nmk$par[7], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[3]), 
                Ec50s     = c(x[4], x[5]), 
                a         = x[6], b = x[7],
                interact  = interact
              )
            }, 
            upper   = upper, 
            lower   = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2:3]), 
            Ec50s  = c(res_IA_nmk$par[4:5]), 
            a      = res_IA_nmk$par[6], 
            b      = res_IA_nmk$par[7], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End DL
    } 
    ## 2.2. Common slopes ----
    else {
      
      ### 2.2.1. IA model ----
      if (interact == "none") {
        
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10)
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat)
        }
        
        if (error_type == "Normal") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3], x[4]), 
                interact  = interact
              )
            }, 
            upper = upper, 
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3], x[4]), 
                interact  = interact
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3], x[4]), 
                interact  = interact
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End IA
      
      ### 2.2.2. SA model ----
      if (interact == "SA") {
        
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10, 20)
        } # the intervals for a and b must be larger beacuse they can compensate each other
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, -20)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat, 0)
        }
        
        if (error_type == "Normal") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], Slopes = c(x[2], x[2]), 
                Ec50s     = c(x[3], x[4]), a = x[5], 
                interact  = interact
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]),
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            a      = res_IA_nmk$par[5], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_Poisson(
                C_mat = C_mat, 
                Response = Response, 
                Max = x[1], 
                Slopes = c(x[2], x[2]), 
                Ec50s = c(x[3], x[4]), 
                a = x[5], 
                interact = interact
              )
            }, 
            upper = upper, 
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s = c(res_IA_nmk$par[3:4]), 
            a = res_IA_nmk$par[5],
            Error = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_NB(
                C_mat = C_mat, 
                Response = Response, 
                Max = x[1], 
                Slopes = c(x[2], x[2]), 
                Ec50s = c(x[3], x[4]), 
                a = x[5], 
                interact = interact
              )
            }, 
            upper = upper, 
            lower = lower,
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s = c(res_IA_nmk$par[3:4]), 
            a = res_IA_nmk$par[5],
            Error = res_IA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End SA
      
      ### 2.2.3. DR model ----
      if (interact == "DR") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10, 20, rep(30, (dim(C_mat)[2] - 1)))
        } # the intervals for a and b must be larger beacuse they can compensate each other
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, -20, rep(-30, (dim(C_mat)[2] - 1)))
        }
        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat, 0, 0)
        }
        
        if (error_type == "Normal") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[2]),
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6:length(x)],
                interact  = interact
              )
            }, 
            upper = upper,
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1],
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            a      = res_IA_nmk$par[5],
            b      = res_IA_nmk$par[6:(6 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6:length(x)], 
                interact  = interact
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            a      = res_IA_nmk$par[5],
            b      = res_IA_nmk$par[6:(6 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response, 
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6:length(x)], 
                interact  = interact
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            a      = res_IA_nmk$par[5],
            b      = res_IA_nmk$par[6:(6 + (dim(C_mat)[2] - 1) - 1)], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      }
      
      ### 2.2.4. DL model ----
      if (interact == "DL") {
        if (missing(upper)) {
          upper <- c(max(Response) * 5, 100, mean_C_mat * 10, 20, 30)
        }
        if (missing(lower)) {
          lower <- c(0, 0, 0, 0, -20, -30)
        }
        if (missing(start)) {
          start <- c(max(Response), 1, mean_C_mat, 0, 0)
        }
        
        if (error_type == "Normal") {
          res_IA_nmk <- nmkb(
            par = start, fn = function(x) {
              IA_complete2_RSS(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]),
                a         = x[5], 
                b         = x[6], 
                interact  = interact
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]),
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            a      = res_IA_nmk$par[5], 
            b      = res_IA_nmk$par[6],
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "Poisson") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_Poisson(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6], 
                interact  = interact
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            a      = res_IA_nmk$par[5],
            b      = res_IA_nmk$par[6], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else if (error_type == "NB") {
          res_IA_nmk <- nmkb(
            par = start, 
            fn = function(x) {
              IA_complete2_NB(
                C_mat     = C_mat, 
                Response  = Response,
                Max       = x[1], 
                Slopes    = c(x[2], x[2]), 
                Ec50s     = c(x[3:4]), 
                a         = x[5], 
                b         = x[6], 
                interact  = interact
              )
            }, 
            upper = upper, 
            lower = lower, 
            control = list(tol = tol_drc)
          )
          Res <- list(
            Max    = res_IA_nmk$par[1], 
            Slopes = c(res_IA_nmk$par[2], res_IA_nmk$par[2]), 
            Ec50s  = c(res_IA_nmk$par[3:4]), 
            a      = res_IA_nmk$par[5],
            b      = res_IA_nmk$par[6], 
            Error  = res_IA_nmk$value
          )
          return(Res)
        } else {
          stop("Misspecification of the error model")
        }
      } # End DL
    } # End common slopes
  } # End Dose-response curves not known
} # End.





f_Jonker_CI_calculations <- function(df_data, N_draws, interact, error_type){
  
  # Pick a number btw 1 and N_obs for each draw
  min_resampling <- 1
  N_obs <- length(df_data$Response)
  max_resampling <- N_obs
  
  # Draws
  draws <- lapply(
    1:N_draws, 
    function(i) {
      n <- sample(min_resampling:max_resampling, 1) 
      lines <- sample(1:N_obs, size = n, replace = TRUE)  
      list(ID = i, N_resampling = n, Lines = lines)
    }
  )
  
  df_draws <- do.call(rbind, lapply(draws, as.data.frame))
  df_draws <- as.data.frame(df_draws)
  df_draws$Lines <- lapply(df_draws$Lines, unlist) 
  
  List_a <- vector("list", N_draws)
  List_b <- vector("list", N_draws)
  
  res <- data.frame(
    a_2.5 = NA, 
    a_97.5 = NA,
    b_2.5 = NA, 
    b_97.5 = NA
    )
  
  i <- 1
  
  for(draw_i in unique(df_draws$ID)){
    print(i)
    
    Lines_i <- unlist(subset(df_draws, ID==draw_i)$Lines)
    
    df_data_i <- df_data |> 
      dplyr::slice(Lines_i)
    
    mod <- tryCatch({
      # Code that can fail
      CA_complete_fit_speed(                         
        C_mat=cbind(                                              
          df_data_i$Dose_EPX,                                  
          df_data_i$Dose_IMD                                  
        ),                                                         
        Response=df_data_i$Response,                         
        interact=interact,                                            
        param=param,
        error_type = error_type,
        multicore = TRUE,                                       
        mc.cores = 4   
      )        
    }, error = function(e) {
      NA
    })
    
    if (!is.na(mod)[1]){
      if (interact == "SA"){
        List_a[[i]] <- mod$a
      }else if (interact %in% c("DL", "DR")){
        List_a[[i]] <- mod$a
        List_b[[i]] <- mod$b
      }else{
        warning("Interaction misspecification")
      }
      
    }else{
      List_a[[i]] <- NA
      List_b[[i]] <- NA
    }
    i <- i+1
  }
  
  if (interact=="SA"){
    res$a_2.5 <- quantile(unlist(List_a), probs = c(0.025, 0.975), na.rm = TRUE)[1]
    res$a_97.5 <- quantile(unlist(List_a), probs = c(0.025, 0.975), na.rm = TRUE)[2]
  }else if(interact %in% c("DL", "DR")){
    res$a_2.5 <- quantile(unlist(List_a), probs = c(0.025, 0.975), na.rm = TRUE)[1]
    res$a_97.5 <- quantile(unlist(List_a), probs = c(0.025, 0.975), na.rm = TRUE)[2]
    res$b_2.5 <- quantile(unlist(List_b), probs = c(0.025, 0.975), na.rm = TRUE)[1]
    res$b_97.5 <- quantile(unlist(List_b), probs = c(0.025, 0.975), na.rm = TRUE)[2]
  }
  
  return(res)
  
}

f_Jonker_CI_calculations_condition <- function(df_data, N_draws, interact, reference = "CA", error_type, max_per_cond = 4){
  
  df_data <- df_data |> 
    mutate(
      condition = paste0(Line, Ratio)
    )
  
  # Nombres d'observations par condition
  cond_levels <- unique(df_data$condition)
  N_obs_cond  <- table(df_data$condition)
  
  List_a <- vector("list", N_draws)
  List_b <- vector("list", N_draws)
  
  res <- data.frame(
    a_2.5 = NA, 
    a_97.5 = NA,
    b_2.5 = NA, 
    b_97.5 = NA
  )

  for(i in seq_len(N_draws)){
    print(i)
    
    # Stratified draws per condition
    n_cond <- sample(1:max_per_cond, 1)
    
    Lines_i <- unlist(lapply(cond_levels, function(cond){
      n_max <- min(N_obs_cond[[cond]], n_cond) 
      sample(
        which(df_data$condition == cond), 
        size = n_max, 
        replace = TRUE
      )
    }))
    
    df_data_i <- dplyr::slice(df_data, Lines_i)
    
    if (reference == "CA") {
      mod <- tryCatch({
        CA_complete_fit_speed(                         
          C_mat=cbind(
            df_data_i$Dose_EPX,
            df_data_i$Dose_IMD
          ),                                                         
          Response=df_data_i$Response,                         
          interact=interact,                                            
          param=param,
          error_type = error_type,
          multicore = TRUE,                                       
          mc.cores = 4   
        )        
      }, error = function(e) NA)
    } else if (reference == "IA") {
      mod <- tryCatch({
        IA_complete_fit_speed(                         
          C_mat=cbind(
            df_data_i$Dose_EPX,
            df_data_i$Dose_IMD
          ),                                                         
          Response=df_data_i$Response,                         
          interact=interact,                                            
          param=param,
          error_type = error_type
        )        
      }, error = function(e) NA)
    }
    
    if (!is.na(mod)[1]){
      if (interact == "SA"){
        List_a[[i]] <- mod$a
      } else if (interact %in% c("DL", "DR")){
        List_a[[i]] <- mod$a
        List_b[[i]] <- mod$b
      } else {
        warning("Interaction misspecification")
      }
    } else {
      List_a[[i]] <- NA
      List_b[[i]] <- NA
    }
  }
  
  if (interact=="SA"){
    res$a_2.5  <- quantile(unlist(List_a), probs = 0.025, na.rm = TRUE)
    res$a_97.5 <- quantile(unlist(List_a), probs = 0.975, na.rm = TRUE)
  } else if (interact %in% c("DL", "DR")){
    res$a_2.5  <- quantile(unlist(List_a), probs = 0.025, na.rm = TRUE)
    res$a_97.5 <- quantile(unlist(List_a), probs = 0.975, na.rm = TRUE)
    res$b_2.5  <- quantile(unlist(List_b), probs = 0.025, na.rm = TRUE)
    res$b_97.5 <- quantile(unlist(List_b), probs = 0.975, na.rm = TRUE)
  }
  
  return(res)
}



