library(here)
source(file = here::here("functions/fun.R"))
f_load_libraries_colors()

# 1. Load data ----
df_mix_coc <- read_excel(here::here("data/Data_Mix.xlsx"), sheet="ReproSpring2025coc") |> 
  mutate(
    Condition = paste0(Ratio, Line),
    w_coc = w_coc_tot/(Nb_cocoons-Nb_cocoons_crushed),
    TU_drc = Dose_IMD/IMD_EC50_drc + Dose_EPX/EPX_EC50_drc,
    
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

Mean_coc_controls <- mean(subset(df_mix_coc, Dose_IMD == 0 & Dose_EPX == 0)$Nb_cocoons)

Shift <- 1e-6

df_mix_coc <- df_mix_coc |> 
  mutate(
    Nb_cocoons = if_else(Nb_cocoons == 0, Shift, Nb_cocoons)
  )

# 2. Entrées du modèle ----

df_design <- read.csv(here::here("data/Design_mixture.csv")) |> 
  filter(EPX == 0)


# df_mix_coc_RatioI <- df_mix_coc |> 
#   mutate(Response = Nb_cocoons) |> 
#   dplyr::select(Dose_EPX, Dose_IMD, Response, Line, Ratio) |> 
#   filter(Ratio == "I")
# 
# df_mix_coc_RatioE <- df_mix_coc_RatioI |> 
#   mutate(Ratio = "E")
# 
# df_mix_coc_int <- rbind(df_mix_coc_RatioI, df_mix_coc_RatioE) |> 
#   mutate(
#     Dose_EPX = if_else(Ratio == "E", Dose_IMD, 0),
#     Dose_IMD = if_else(Ratio == "E", 0, Dose_EPX)
#     )


df_mix_coc_int <- df_mix_coc |> 
  mutate(Response = Nb_cocoons) |> 
  dplyr::select(Dose_EPX, Dose_IMD, Response, Line, Ratio) |> 
  as.data.frame() |> 
  filter(!(Dose_IMD == 0 & Dose_EPX == 0))

df_mix_drc <- df_mix_coc_int |> 
  filter(Ratio %in% c("I", "E")) |> 
  mutate(
    Dose = if_else(Ratio == "I", Dose_IMD, Dose_EPX)
  )


drc.mix<- drm(
  Response~Dose, Ratio, 
  data = df_mix_drc,
  type="Poisson",
  fct=LL.4(
    names = c("slope", "Ymin", "Ymax", "EC50"),
    fixed = c(NA, 0, NA, NA)
  ), 
  # Reproduction is expected to reach 0 for large concentrations
  pmodels=data.frame(1, 1, Ratio)
)

summary(drc.mix)


C_mat <- cbind(                                                    
  df_mix_drc$Dose_EPX,                                      
  df_mix_drc$Dose_IMD                                       
)

mean_C_mat <- c(mean(df_mix_drc$Dose_EPX), mean(df_mix_drc$Dose_IMD))

Response <- df_mix_drc$Response

multicore = TRUE

mc.cores <- 4

upper=c(max(Response)*2, 100, mean_C_mat*10)
lower=c(0,0,0,0)
start=c(max(Response),2,mean_C_mat)

res_CA_nmk_poisson<-nmkb(
  par=start, 
  fn=function(x){
    CA_complete2_Poisson(
      C_mat=C_mat, 
      Response=Response, 
      Max=x[1], 
      Slopes=c(x[2], x[2]), 
      Ec50s=c(x[3], x[4]), 
      interact=interact, 
      multicore=multicore, 
      mc.cores=mc.cores
      )
    }, 
  upper=upper, 
  lower=lower, 
  control=list(tol=10^-8)
  )

res_CA_nmk_normal<-nmkb(
  par=start, 
  fn=function(x){
    CA_complete2_RSS(
      C_mat=C_mat, 
      Response=Response, 
      Max=x[1], 
      Slopes=c(x[2], x[2]), 
      Ec50s=c(x[3], x[4]), 
      interact=interact, 
      multicore=multicore, 
      mc.cores=mc.cores
    )
  }, 
  upper=upper, 
  lower=lower, 
  control=list(tol=10^-8)
)

res_CA_nmk_poisson
res_CA_nmk_normal

# Do not work with the controls in the data frame !!!!!!
CA_fit<-CA_complete2_fit_speed(                                    
  C_mat=cbind(                                                    
    df_mix_coc_int$Dose_EPX,                                      
    df_mix_coc_int$Dose_IMD                                       
  ),                                                            
  Response=df_mix_coc_int$Response,                               
  interact="none",                                                  
  # param=param,
  error_type = "Poisson",
  identical_slopes = TRUE,
  multicore = TRUE
)
CA_fit



C_mat <- cbind(                                                    
  df_mix_coc_int$Dose_EPX,                                      
  df_mix_coc_int$Dose_IMD                                       
)  

Response <- df_mix_coc_int$Response                               
interact <- "none"                                                  
# param=param,
error_type <- "Poisson"
identical_slopes <- TRUE

param <- NULL



if (is.null(param)){
  
  require(dfoptim)
  
  if (any(C_mat==0)){
    mean_C_mat<-exp(apply(log(C_mat[-which(C_mat==0, arr.ind=TRUE)[,1],]),2,mean))
  }else{
    mean_C_mat<-exp(apply(log(C_mat),2,mean))
  }
  
  mean_C_mat <- 0.5
  
  if (identical_slopes) {
    
    if (interact=="none"){
      
      if (missing(upper))
        upper=c(max(Response)*5, 100, mean_C_mat*10)
      
      if (missing(lower))
        lower=c(0,0,0,0)
      
      if (missing(start))
        start=c(max(Response),1,mean_C_mat)
      
      # Fit
      if (error_type == "Normal"){
        res_CA_nmk<-nmkb(par=start, 
                         fn=function(x){CA_complete2_RSS(C_mat=C_mat, 
                                                         Response=Response, 
                                                         Max=x[1], 
                                                         Slopes=c(x[2], x[2]), 
                                                         Ec50s=c(x[3], x[4]), 
                                                         interact=interact, 
                                                         multicore=multicore, 
                                                         mc.cores=mc.cores)}, 
                         upper=upper, 
                         lower=lower, 
                         control=list(tol=10^-6))
      }else if(error_type == "Poisson"){
        res_CA_nmk<-nmkb(par=start, 
                         fn=function(x){CA_complete2_Poisson(C_mat=C_mat, 
                                                             Response=Response, 
                                                             Max=x[1], 
                                                             Slopes=c(x[2], x[2]), 
                                                             Ec50s=c(x[3], x[4]), 
                                                             interact=interact, 
                                                             multicore=multicore, 
                                                             mc.cores=mc.cores)}, 
                         upper=upper, 
                         lower=lower, 
                         control=list(tol=10^-6))
      }else{
        stop("Misspecification of the error model type")
      }
      
      return(list(Max=res_CA_nmk$par[1], 
                  Slopes=c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
                  Ec50s=c(res_CA_nmk$par[3:4]), 
                  RSS=res_CA_nmk$value))
    }
  }
}

CA_complete2 <- function(C_mat, Max, Slopes, Ec50s, a=0, b=0, interact="none", multicore=FALSE, mc.cores=4){
  param<-data.frame(Slopes, Ec50s)
  #try(cat(Slopes))
  #try(cat(Ec50s))
  if (interact=="none"){
    a<-0
    b<-0
    f<-function(Y, x=x, Max=Max, param=param, a=a, b=b){
      ecs <- param$Ec50s * ((Max - Y)/Y) ^ (1 / param$Slopes)  # C1 corresponding to Y 
      G<-(sum(x / ecs))-1
      return(abs(G))
    }
  }
  else if (interact=="SA"){
    b<-0
    f<-function(Y, x=x, Max=Max, param=param, a=a, b=b){
      ecs <- param$Ec50s * ((Max - Y)/Y) ^ (1 / param$Slopes)  # C1 corresponding to Y 
      G<-(sum(x / ecs))-exp(a*prod(x/param$Ec50s/(sum(x/param$Ec50s))))
      return(abs(G))
    }
  }
  else if (interact=="DR"){
    if (length(b)!=(dim(C_mat)[2]-1))
      stop("length of b should be equal to the number of chemicals -1.")
    f<-function(Y, x=x, Max=Max, param=param, a=a, b=b){
      ecs <- param$Ec50s * ((Max - Y)/Y) ^ (1 / param$Slopes)  # Cs corresponding to Y 
      G<-(sum(x / ecs))-exp((a+b%*%((x/param$Ec50s/(sum(x/param$Ec50s))))[1:(dim(C_mat)[2]-1)])*prod(x/param$Ec50s/(sum(x/param$Ec50s))))
      return(abs(G))
    }
  }
  else if (interact=="DR2"){
    if (length(b)!=(dim(C_mat)[2]-1))
      stop("length of b should be equal to the number of chemicals -1.")
    f<-function(Y, x=x, Max=Max, param=param, a=a, b=b){
      ecs <- param$Ec50s * ((Max - Y)/Y) ^ (1 / param$Slopes)  # Cs corresponding to Y 
      G<-(sum(x / ecs))-exp((a+b%*%sin((x/param$Ec50s/(sum(x/param$Ec50s)))*2*pi)[1:(dim(C_mat)[2]-1)])*prod(x/param$Ec50s/(sum(x/param$Ec50s))))
      return(abs(G))
    }
  }
  else if (interact=="DL"){
    f<-function(Y, x=x, Max=Max, param=param, a=a, b=b){
      if (length(b)!=(dim(C_mat)[2]-1))
        stop("length of b should be equal to the number of chemicals -1.")
      ecs <- param$Ec50s * ((Max - Y)/Y) ^ (1 / param$Slopes)  # C1 corresponding to Y 
      G<-(sum(x / ecs))-exp(a*(1-b*(sum(x/param$Ec50s)))*prod(x/param$Ec50s/(sum(x/param$Ec50s))))
      return(abs(G))
    }
  }
  else 
    stop("please specify interaction model")
  Y_f<-function(x, Max=Max, param=param, a=a, b=b) {
    if (all(x==0)){          #no chemical
      if  (all(param$Slopes > 0))
        return(Max)
      if  (all(param$Slopes < 0))
        return(0)
    }
    #Y<-optimize(f=function(Y) f(Y=Y, x=x, Max=Max, param=param, a=a, b=b), interval=c(0,Max*1.2))$minimum
    Y<-optimize(f=function(Y) f(Y=Y, x=x, Max=Max, param=param, a=a, b=b), interval=c(0,Max))$minimum
    
    return(Y)
  }
  if (multicore){
    C_mat_list<-split(unique(C_mat), cut(1:(dim(unique(C_mat))[1]), breaks=dim(unique(C_mat))[1]))
    res_CA_unique<-unlist(mclapply(C_mat_list, function(x) Y_f(x,Max=Max, param=param, a=a, b=b), mc.cores=mc.cores))
    res_CA<-rep(NA, dim(unique(C_mat))[1])
    for (i in 1:(dim(unique(C_mat))[1]))
      res_CA[which((C_mat[,1]==unique(C_mat)[i,1])&(C_mat[,2]==unique(C_mat)[i,2]))]<-res_CA_unique[i]
  }
  else {
    res_CA_unique<-apply(unique(C_mat),1, function(x) Y_f(x, Max=Max, param=param, a=a, b=b))
    res_CA<-rep(NA, dim(unique(C_mat))[1])
    for (i in 1:(dim(unique(C_mat))[1]))
      res_CA[which((C_mat[,1]==unique(C_mat)[i,1])&(C_mat[,2]==unique(C_mat)[i,2]))]<-res_CA_unique[i]
  }
  if (any(is.na(res_CA)))
    cat(param, a, b)
  return(res_CA)
} 

CA_complete2_fit_speed<-function(C_mat, Response, param=NULL, upper=NULL, lower=NULL, start=NULL, interact="none", identical_slopes=FALSE, error_type = "Normal", iter=500, multicore=FALSE, mc.cores=4){
  
  if (is.null(param)){
    
    require(dfoptim)
    
    if (any(C_mat==0))
      mean_C_mat<-exp(apply(log(C_mat[-which(C_mat==0, arr.ind=TRUE)[,1],]),2,mean))
    
    else
      mean_C_mat<-exp(apply(log(C_mat),2,mean))
    
    if (identical_slopes) {
      
      if (interact=="none"){
        
        if (missing(upper))
          upper=c(max(Response)*5, 0, mean_C_mat*10)
        
        if (missing(lower))
          lower=c(0,-100,0,0)
        
        if (missing(start))
          start=c(max(Response),-1,mean_C_mat)
        
        # Fit
        if (error_type == "Normal"){
          res_CA_nmk<-nmkb(par=start, 
                           fn=function(x){CA_complete2_RSS(C_mat=C_mat, 
                                                           Response=Response, 
                                                           Max=x[1], 
                                                           Slopes=c(x[2], x[2]), 
                                                           Ec50s=c(x[3], x[4]), 
                                                           interact=interact, 
                                                           multicore=multicore, 
                                                           mc.cores=mc.cores)}, 
                           upper=upper, 
                           lower=lower, 
                           control=list(tol=10^-6))
        }else if(error_type == "Poisson"){
          res_CA_nmk<-nmkb(par=start, 
                           fn=function(x){CA_complete2_Poisson(C_mat=C_mat, 
                                                               Response=Response, 
                                                               Max=x[1], 
                                                               Slopes=c(x[2], x[2]), 
                                                               Ec50s=c(x[3], x[4]), 
                                                               interact=interact, 
                                                               multicore=multicore, 
                                                               mc.cores=mc.cores)}, 
                           upper=upper, 
                           lower=lower, 
                           control=list(tol=10^-6))
        }else{
          stop("Misspecification of the error model type")
        }
        
        return(list(Max=res_CA_nmk$par[1], 
                    Slopes=c(res_CA_nmk$par[2], res_CA_nmk$par[2]), 
                    Ec50s=c(res_CA_nmk$par[3:4]), 
                    RSS=res_CA_nmk$value))
      }
    }
  }
}
