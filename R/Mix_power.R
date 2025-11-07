library(tidyverse)
library(here)

source(file = "fun.R")

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
param <- data.frame(Slopes=Slopes, Max=Max, Ec50s=Ec50s)


f_power_design <- function(Nb_rep, l_interaction_a, Nb_pop, param){
  
  res <- data.frame()
  
  for (Nb_rep_j in Nb_rep){
    for (interaction_a_i in l_interaction_a){
      
     # print(interaction_a_i)
      
      a_i <- interaction_a_i
      
      # Calculation of the theorical response
      res_CA_theorical_design_i <- CA_complete2(
        C_mat, 
        Max, Slopes, Ec50s, 
        interact="SA", 
        a=a_i
      )
      
      # We add the corresponding response to our design
      df_dose_rep_i <- df_design%>%
        mutate(rep_Y=res_CA_theorical_design_i)
      
      df_dose_rep_i <- df_dose_rep_i[,-c(1, 4)]
      
      for(k in 1:Nb_pop){
        
        set.seed(121212+k)
        ID_pop <- k
        # Population simulation 
        df_sim_data_k <- data.frame(
          Y=c(),
          EPX=c(),
          IMD=c()
        )
        
        for(i in 1:length(df_dose_rep_i$EPX)){
          
          Y_th <- df_dose_rep_i$rep_Y[i]
          Y_sim <- rpois(Nb_rep_j, Y_th) # <1>
          
          EPX <- rep(df_dose_rep_i$EPX[i], Nb_rep_j)
          IMD <- rep(df_dose_rep_i$IMD[i], Nb_rep_j)
          
          df_sim_data_j <- data.frame(
            Y_th=rep(Y_th, Nb_rep_j),
            Y_sim=Y_sim,
            EPX=EPX,
            IMD=IMD
          )
          df_sim_data_k <- rbind(df_sim_data_k, df_sim_data_j)
        }
        # Fit of the model 
        res_CA_SA_fit_k<-CA_complete2_fit_speed(
          C_mat=cbind(                      
            df_sim_data_k$EPX,                
            df_sim_data_k$IMD                 
          ),                              
          Response=df_sim_data_k$Y_sim,       
          interact="SA",                      
          param=param,
          error_type = "Poisson"
        )
        res_a_k <- res_CA_SA_fit_k$a
        
        res_k <- data.frame(
          Simulated_interaction = a_i,
          Estimated_interaction = res_a_k,
          ID_pop = ID_pop,
          Nb_rep = Nb_rep_j
        )
        res <- rbind(res, res_k)
      }
      
      
      #print("Calculations done for given interaction")
      
    }
    #print("Calculations done for given Nb_rep")
  }
  
  return(res)
}

Nb_rep <- c(4)
tmp <- c(10, 7, 5, 2, 1)
l_interaction_a <- c(tmp, -rev(tmp))
Nb_pop <- 10

df_boot_res <- f_power_design(Nb_rep, l_interaction_a, Nb_pop, param)
save(df_boot_res,file = here::here("mod/MIX/Res_boot_design_new.RData"))


