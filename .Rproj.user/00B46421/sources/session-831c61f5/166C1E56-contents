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
      } else {
        stop("Misspecification of the error model")
      }
    } # End SA model
  } # End dose-response curves known

  # 2. If dose-response curves known ----
  if (is.null(param)) {
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
            upper   = upper,
            lower   = lower,
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
  res_CA_Poisson <- -sum(Response * log(res_CA) - res_CA - lfactorial(Response)) # - Loglikelihood

  return(res_CA_Poisson)
}

CA_complete_fit_speed <- function(C_mat, Response, param = NULL, upper = NULL, lower = NULL, start = NULL, interact = "none", identical_slopes = FALSE, error_type = "Normal", iter = 500, multicore = FALSE, mc.cores = 4) {
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
          upper <- c(max(Response) * 5, 0, 0, mean_C_mat * 10)
        }
        if (missing(lower)) {
          lower <- c(0, -100, -100, 0, 0)
        }
        if (missing(start)) {
          start <- c(max(Response), -1, -1, mean_C_mat)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
          upper <- c(max(Response) * 5, 0, 0, mean_C_mat * 10, 20) # the intervals for a and b must be larger beacuse they can compensate each other
        }
        if (missing(lower)) {
          lower <- c(0, -100, -100, 0, 0, -20)
        }
        if (missing(start)) {
          start <- c(max(Response), -1, -1, mean_C_mat, 0)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
          upper <- c(max(Response) * 5, 0, 0, mean_C_mat * 10, 20, rep(30, (dim(C_mat)[2] - 1)))
        } # the intervals for a and b must be larger beacuse they can compensate each other
        if (missing(lower)) {
          lower <- c(0, -100, -100, 0, 0, -20, rep(-30, (dim(C_mat)[2] - 1)))
        }
        if (missing(start)) {
          start <- c(max(Response), -1, -1, mean_C_mat, 0, rep(0, (dim(C_mat)[2] - 1)))
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
          upper <- c(max(Response) * 5, 0, 0, mean_C_mat * 10, 20, 30)
        }
        if (missing(lower)) {
          lower <- c(0, -100, -100, 0, 0, -20, -30)
        }
        if (missing(start)) {
          start <- c(max(Response), -1, -1, mean_C_mat, 0, 0)
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
            control = list(tol = 10^-6)
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
              control = list(tol = 10^-6)
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
          upper <- c(max(Response) * 5, 0, mean_C_mat * 10)
        }
        if (missing(lower)) {
          lower <- c(0, -100, 0, 0)
        }
        if (missing(start)) {
          start <- c(max(Response), -1, mean_C_mat)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
          upper <- c(max(Response) * 5, 0, mean_C_mat * 10, 20)
        } # the intervals for a and b must be larger beacuse they can compensate each other
        if (missing(lower)) {
          lower <- c(0, -100, 0, 0, -20)
        }
        if (missing(start)) {
          start <- c(max(Response), -1, mean_C_mat, 0)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
          upper <- c(max(Response) * 5, 0, mean_C_mat * 10, 20, rep(30, (dim(C_mat)[2] - 1)))
        } # the intervals for a and b must be larger beacuse they can compensate each other
        if (missing(lower)) {
          lower <- c(0, -100, 0, 0, -20, rep(-30, (dim(C_mat)[2] - 1)))
        }
        if (missing(start)) {
          start <- c(max(Response), -1, mean_C_mat, 0, 0)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
          upper <- c(max(Response) * 5, 0, mean_C_mat * 10, 20, 30)
        }
        if (missing(lower)) {
          lower <- c(0, -100, 0, 0, -20, -30)
        }
        if (missing(start)) {
          start <- c(max(Response), -1, mean_C_mat, 0, 0)
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
            control = list(tol = 10^-6)
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
            control = list(tol = 10^-6)
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
