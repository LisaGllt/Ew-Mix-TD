
# Fonctions - Mixture - Jonker interaction models ----

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

# Application du modèle ----

## Simulation d'un jeu de données ----

Nb_rep <- 4
# Simulated interaction
interaction_a <- -1

df_design <- read.csv(here::here("data/Design_mixture.csv"))
df_drc_param <- read.csv(here::here("data/DRC_parameters_cocoons_really_used.csv"))

slope <- df_drc_param$coef.drc.4.[1]
Y_max <- df_drc_param$coef.drc.4.[2]
EPX_EC50 <- df_drc_param$coef.drc.4.[3]
IMD_EC50 <- df_drc_param$coef.drc.4.[4]
y_min <- 0

C_mat = cbind(df_design$EPX,df_design$IMD)
Max = Y_max
Slopes = c(slope, slope)
Ec50s = c(EPX_EC50, IMD_EC50)

# Calculation of effects for the design concentrations 
res_CA_theorical_design <- CA_complete2(
  C_mat, 
  Max, Slopes, Ec50s, 
  interact="SA", 
  a=interaction_a
)

# We add the corresponding theorical response to our design
df_dose_rep <- df_design%>%
  mutate(rep_Y=res_CA_theorical_design)

df_dose_rep <- df_dose_rep[,-c(1, 4)]

df_sim_data <- data.frame(
  Y=c(),
  EPX=c(),
  IMD=c()
)

set.seed(12)
for(i in 1:length(df_dose_rep$EPX)){
  
  Y_th <- df_dose_rep$rep_Y[i] # Theorical response for a given condition
  Y_sim <- rpois(Nb_rep, Y_th) # <1> # Simulated responses for a given condition
  
  EPX <- rep(df_dose_rep$EPX[i], Nb_rep)
  IMD <- rep(df_dose_rep$IMD[i], Nb_rep)
  
  df_sim_data_i <- data.frame(
    Y_th=rep(Y_th, Nb_rep),       
    Y_sim=Y_sim,                  
    EPX=EPX,
    IMD=IMD
  )
  df_sim_data <- rbind(df_sim_data, df_sim_data_i)
}

## Fit du modèle sur le jeu de données ----

param <- data.frame(Slopes=Slopes, Max=Max, Ec50s=Ec50s)

df_data_fit<-CA_complete2_fit_speed(
  C_mat=cbind(                     
    df_sim_data$EPX,               
    df_sim_data$IMD                 
  ),                            
  Response=df_sim_data$Y_sim,        
  interact="SA",                    
  param=param,                    
  multicore = TRUE               
)                                 

# Log-likelihood of the model
LL_data_fit_SA <- df_data_fit$Error

res_a <- round(df_data_fit$a, 3)                          
ratio_interaction_a=round(df_data_fit$a/interaction_a, 3)

df_res <- data.frame(`Simulated interaction parameter`=interaction_a,
                     `Estimated interaction parameter`=round(df_data_fit$a,3),
                     Ratio=ratio_interaction_a,
                     LL = LL_data_fit_SA)
df_res |> datatable(options = list(dom = 't'), class="hover") 

## Représentation graphique ----

# Graphic preparation : Construction of the surface mesh

min_dose_EPX <- min(subset(df_design, !EPX==0)$EPX)
max_dose_EPX <- max(subset(df_design, !EPX==0)$EPX)

min_dose_IMD <- min(subset(df_design, !IMD==0)$IMD)
max_dose_IMD <- max(subset(df_design, !IMD==0)$IMD)

by_EPX <- (log(max_dose_EPX)-log(min_dose_EPX))/100
by_IMD <- (log(max_dose_IMD)-log(min_dose_IMD))/100

# Grid creation
grid_C1 <- exp(seq(log(min_dose_EPX), log(max_dose_EPX), by = by_EPX))
grid_C2 <- exp(seq(log(min_dose_IMD), log(max_dose_IMD), by = by_IMD))

grid <- expand.grid(
  x = grid_C1, 
  y = grid_C2
)

df_data_solved_surf <- CA_complete2(
  C_mat=grid,                    
  Max=Max,                        
  Slopes=Slopes,                
  Ec50s=Ec50s,                      
  a=df_data_fit$a, interact="SA"   
)                                  

df_data_solved_surf_ordered <- df_data_solved_surf[order(grid[,2])]

z_SA_Est <- matrix(
  df_data_solved_surf_ordered,
  nrow=length(grid_C1),
  ncol=length(grid_C2),
  byrow=FALSE
)

grid_EstSurface_SA <- expand.grid(x = grid_C1, y = grid_C2)
grid_EstSurface_SA$z <- as.vector(z_SA_Est)


target_effect <- 50
bin_contour <- 1

p <- ggplot() +
  geom_contour(
    data = grid_CA,
    aes(x = x, y = y, z = z),
    color = "black",
    linewidth = 0.8,
    linetype = "dotted",
    binwidth = bin_contour
  ) +
  geom_contour(
    data = grid_ThSurface_SA,
    aes(x = x, y = y, z = z),
    color = "blue4", 
    linewidth = 1,
    binwidth = bin_contour
  ) +
  metR::geom_text_contour(
    data = grid_ThSurface_SA,
    aes(x = x, y = y, z = z), 
    stroke = 0.2,
    size = 4,
    color = "blue4",
    skip = 1,
    label.placer = label_placer_flattest(),
    binwidth = bin_contour
  ) +
  labs(
    x = expression("[EPX] (mg·kg"^-1*")"),
    y = expression("[IMD] (mg·kg"^-1*")")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

p
