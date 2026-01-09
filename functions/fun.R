# 1. General setup ----

# LOAD ALL REQUIRED LIBRAIRIES
f_load_libraries <- function() {
  
  # 📦 Data manipulation & I/O
  library(tidyverse)     # Data manipulation and visualization
  library(data.table)    # Fast data manipulation
  library(here)          # File path management
  library(readxl)        # Read Excel files
  library(reshape2)      # Data reshaping
  library(DT)            # Interactive tables
  library(knitr)         # Dynamic reports (RMarkdown / Quarto)
  
  # 🎨 Visualization & graphics
  library(ggplot2)       # Data visualization
  library(ggthemes)      # ggplot2 themes
  library(ggdist)        # Distributions and uncertainty
  library(ggsci)         # Scientific color palettes
  library(viridis)       # Perceptually uniform color palette
  library(wesanderson)   # Artistic color palettes
  library(RColorBrewer)  # Predefined color palettes
  library(nord)          # Nord-inspired color palettes
  library(colorspace)    # Color spaces
  library(scales)        # ggplot2 scales
  
  library(plotly)        # Interactive graphics
  library(ggiraph)       # Interactive ggplot2 graphics
  library(ggrepel)       # Non-overlapping text labels
  library(patchwork)     # Combine multiple ggplots
  library(gridExtra)     # Grid-based plot layouts
  library(grid)          # Low-level graphical tools
  library(ggbreak)       # Broken axes in ggplot2
  library(ggtext)        # Advanced text formatting in ggplot2
  library(ggh4x)         # Advanced ggplot2 extensions
  library(metR)          # Meteorological and geoscience tools
  library(rnaturalearth) # Natural Earth spatial data
  library(sf)            # Spatial data handling
  
  # 📋 Tables & reporting
  library(kableExtra)    # Enhanced tables
  library(flextable)     # Word / PowerPoint tables
  library(gt)            # Table rendering
  library(processx)      # Process management
  
  # 📊 Bayesian modeling & inference
  library(brms)          # Bayesian modeling with Stan
  library(rstan)         # R interface to Stan
  library(cmdstanr)      # CmdStan interface
  library(tidybayes)     # Manipulation of Bayesian outputs
  library(ggmcmc)        # MCMC diagnostics
  library(rethinking)    # Advanced Bayesian modeling
  library(priorsense)    # Prior sensitivity analysis
  library(coda)          # MCMC analysis tools
  library(bayesnec)      # Bayesian NEC models
  library(truncnorm)     # Truncated normal distributions
  
  # 🔬 Regression & hypothesis testing
  library(car)           # Regression and statistical tests
  library(nlstools)      # Nonlinear model tools
  library(lsmeans)       # Post-hoc comparisons
  library(ggpubr)        # Publication-ready plots
  library(marginaleffects) # Marginal effects
  library(brglm2)        # Bias-reduced generalized linear models
  library(multcomp)      # Multiple comparisons
  
  # ⚙️ Computation & numerical tools
  library(parallel)      # Parallel computing
  library(deSolve)       # Differential equations
  library(tmvtnorm)      # Truncated multivariate normal
  library(fdrtool)       # False discovery rate
  library(drc)           # Dose–response modeling
  
  # 🧪 Model evaluation & workflow
  library(easystats)     # Statistics and modeling ecosystem
  library(performance)  # Model diagnostics
  library(modelsummary) # Model summaries
  library(plan)         # Task scheduling
  
  # 🖋️ Math, LaTeX & fonts
  library(latex2exp)     # LaTeX expressions in ggplot2
  library(extrafont)     # Font management
}

# LOAD ALL USED AESTHETICS
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

# LOAD LIBRAIRIES AND AESTHETICS
f_load_libraries_colors <- function() {
  f_load_libraries()
  f_load_colors()
}

# 2. Load experimental data ----

# LOAD REPRODUCTION DATA FROM EXPERIMENT A
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

# LOAD HATCHLING DATA FROM EXPERIMENT A
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

# LOAD GROWTH DATA FROM EXPERIMENT A
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
    df_expA_growth_tot = df_growth,
    df_expA_growth = df_growth_alive
  )
  
  return(df_res)
}

# LOAD DATA FROM EXPERIMENT B
f_read_data_expB <- function() {
  
  Val_Ctrl <- 1e-4
  Label_Ctrl <- "0 (Ctrl)"
  
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
  
  df_mix_ind_mean <- df_mix_ind |> 
    aggregate(w ~ Condition, FUN = mean)
  
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
    mutate(
      TU_drc = Dose_IMD / IMD_EC50_drc + Dose_EPX/EPX_EC50_drc,
      Dose_IMD_plot = case_when(
        Dose_IMD == 0 ~ Val_Ctrl,
        !(Dose_IMD == 0) ~ Dose_IMD
      ),
      Dose_EPX_plot = case_when(
        Dose_EPX == 0 ~ Val_Ctrl,
        !(Dose_EPX == 0) ~ Dose_EPX
      )
    )
  
  df_mix_cocsize_cosm <- df_mix_cocsize |> 
    aggregate(v_coc ~ Ratio + Line + Nb_rep + TU_drc, FUN = mean) |> 
    mutate(Condition = paste0(Ratio, Line))
  
  df_mix_coc <- df_mix_coc |> 
    full_join(df_mix_cocsize_cosm, by = c("Condition", "Nb_rep", "Ratio", "Line", "TU_drc")) |> 
    mutate(
      Id_Cond = w_coc/v_coc
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
  
  df_mix_coc_single_mean <- df_mix_coc_single |> 
    aggregate(Nb_cocoons ~ Dose+Molecule+Dose_plot, mean)
  
  Date_D28_LotD <- dmy("23/04/2025")
  Date_D28_LotE <- dmy("29/04/2025")
  
  df_mix_hatch_wide_D <- read_excel(here::here("data/Data_Mix.xlsx"), sheet="ReproSpring2025hatchlingsLotD") |> 
    dplyr::select(-Reste) |> 
    mutate(across(`23/04/2025`:`20/06/2025`, ~replace_na(.x, 0))) |> 
    mutate(
      ID_rep = paste0(Ratio, Line, Nb_rep),
      Condition = paste0(Ratio,Line)
    )
  
  df_mix_hatch_wide_E <- read_excel(here::here("data/Data_Mix.xlsx"), sheet="ReproSpring2025hatchlingsLotE") |> 
    dplyr::select(-Reste)|> 
    mutate(across(`29/04/2025`:`20/06/2025`, ~replace_na(.x, 0))) |> 
    mutate(
      ID_rep = paste0(Ratio, Line, Nb_rep),
      Condition = paste0(Ratio,Line)
    )
  
  df_mix_hatch_long_D <- df_mix_hatch_wide_D %>%
    pivot_longer(
      cols = matches("\\d{2}/\\d{2}/\\d{4}"), 
      names_to = "Date",
      values_to = "Nb_hatch"
    ) %>%
    mutate(Date = dmy(Date)) 
  
  df_mix_hatch_long_E <- df_mix_hatch_wide_E %>%
    pivot_longer(
      cols = matches("\\d{2}/\\d{2}/\\d{4}"), 
      names_to = "Date",
      values_to = "Nb_hatch"
    ) %>%
    mutate(Date = dmy(Date)) 
  
  df_mix_hatch_long <- rbind(df_mix_hatch_long_D, df_mix_hatch_long_E)
  
  df_mix_hatch_long_D_rep <- df_mix_hatch_long_D %>%
    arrange(ID_rep, Date) %>%
    group_by(ID_rep) %>%
    mutate(Cumul_hatch = cumsum(Nb_hatch)) %>%
    ungroup()
  
  df_mix_hatch_long_E_rep <- df_mix_hatch_long_E %>%
    arrange(ID_rep, Date) %>%
    group_by(ID_rep) %>%
    mutate(Cumul_hatch = cumsum(Nb_hatch)) %>%
    ungroup()
  
  df_mix_hatch <- rbind(df_mix_hatch_long_D_rep, df_mix_hatch_long_E_rep) |> 
    mutate(
      t = case_when(
        Lot == "D" ~ as.numeric(Date - Date_D28_LotD),
        Lot == "E" ~ as.numeric(Date - Date_D28_LotE)
      ),
      Cumul_per_hatch = Cumul_hatch/Tot_coc
    )
  
  tot_coc_per_condition <- df_mix_coc %>%
    group_by(Condition) %>%
    summarise(Tot_coc_cond = sum(Nb_cocoons, na.rm = TRUE), .groups = "drop")
  
  df_mix_hatch_mean <- df_mix_hatch_long |> 
    mutate(
      t = case_when(
        Lot == "D" ~ as.numeric(Date - Date_D28_LotD),
        Lot == "E" ~ as.numeric(Date - Date_D28_LotE)
      )
    ) |> 
    arrange(Condition, t, Ratio, Line) |> 
    group_by(Condition, Date, t, Ratio, Line) |> 
    summarise(
      Nb_hatch = sum(Nb_hatch, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(Condition, t, Ratio, Line) |> 
    group_by(Condition) |> 
    mutate(Cumul_hatch = cumsum(Nb_hatch)) |> 
    ungroup() |> 
    left_join(tot_coc_per_condition, by = "Condition") |> 
    mutate(Cumul_perc_hatch = Cumul_hatch / Tot_coc_cond * 100)
  
  
  res <- list(
    df_expB_growth = df_mix_ind,
    df_expB_growth_mean = df_mix_ind_mean,
    df_expB_repro = df_mix_coc,
    df_expB_repro_mean = df_mix_coc_mean,
    df_expB_repro_single = df_mix_coc_single,
    df_expB_repro_single_mean = df_mix_coc_single_mean,
    df_expB_cocsize = df_mix_cocsize,
    df_expB_hatch = df_mix_hatch,
    df_expB_hatch_mean = df_mix_hatch_mean
  )
  
  return(res)
}

# PERFORM A LOG-LIKELIHOOD RATIO TEST BETWEEN A SIMPLE MODEL AND A MORE COMPLEX MODEL
f_LLRatio_test <- function(LL_model, LL_anova, k_model, k_anova, signif) {
  
  # Likelihood ratio test
  # H0 : model is acceptable 
  # H1 : ANOVA model has a better fit
  D <- -2 * (LL_model - LL_anova)
  df_test <- k_anova - k_model
  
  # Vérification
  char_LL <- paste0("LL_anova:", round(LL_anova, signif), " | LL_model:", round(LL_model, signif))
  char_k <- paste0("k_anova:", k_anova, " | k_model:", k_model, " | df:", df_test)
  
  p_val <- pchisq(D, df = df_test, lower.tail = FALSE)
  char_test <- paste0("Likelihood-ratio test: D = ", round(D, signif), ", df = ", df_test, ", p = ", sprintf("%.2e", p_val), " -> Signif ? ", (p_val < 0.05), "\n")
  
  # If p-val < 0.05, we cannot reject H1 = Anova model 
  
  res <- cat(char_LL, "\n", char_k, "\n", char_test)
  
  return(res)
}

# 3. Dose-response curves ----

# RETRIEVE CI OF THE DRC MODEL AND MAKE A TABLE OF ESTIMATED PARAMETERS AND THEIR CI
f_CI_drc <- function(Dose_min, Dose_max, Molecules, drc.model, signif_param) {
  
  Dose_x <- expand.grid(
    exp(
      seq(
        log(Dose_min), log(Dose_max),
        by = (log(Dose_max) - log(Dose_min)) / 100
      )
    )
  )

  df_CI <- NULL

  for (Molecule_i in Molecules) {
    df_CI_i <- data.frame(
      Dose = Dose_x$Var1,
      Molecule = c(rep(Molecule_i, length(Dose_x$Var1)))
    )

    df_CI <- rbind(df_CI, df_CI_i)
  }

  Predictions <- predict(drc.model, newdata = df_CI, interval = "confidence")
  df_CI$p <- Predictions[, 1]
  df_CI$pmin <- Predictions[, 2]
  df_CI$pmax <- Predictions[, 3]

  est <- coef(drc.model)
  ci <- confint(drc.model)
  params_tab <- data.frame(
    Estimate  = unname(est),
    CI_lower  = ci[, 1],
    CI_upper  = ci[, 2]
  ) |>
    rownames_to_column("Parameter") |>
    mutate(
      Value_CI = paste0(
        signif(Estimate, signif_param), " [",
        signif(CI_lower, signif_param), ", ",
        signif(CI_upper, signif_param), "]"
      ),
      Parameter = gsub(":\\(Intercept\\)", "", Parameter)
    ) |>
    dplyr::select(Parameter, Value_CI)

  res <- list(
    df_CI = df_CI,
    df_param = params_tab
  )

  return(res)
}

# 4. Mixture - Jonker interaction models ----

# Functions parameters :

# C_mat : Data frame with the doseA and doseB combinations to describe the surface (doseA and doseB are columns) 
# Response : Response vector corresponding to the response at each surface point from C_mat. 
# Max : Maximum response (numeric)
# Slopes : vector of length 2 of the dose-response curves slopes estimated for each substance
# Ec50s : vector of length 2 of the EC50 estimated for each substance
# a : first interaction parameter (numeric)
# b : second interaction parameter (numeric)
# interact : (character chain)
#   - "none" for CA or IA model
#   - "SA" for CASA or IASA model
#   - "DL" for CADL or IADL model
#   - "DR" for CADR or IADR model
# multicore : use of multivore (bolean)
# mc.cores : number of cores to use (numeric)
# param : param <- data.frame(Slopes = Slopes, Max = Max, Ec50s = Ec50s)
# upper : vector of upper limits for param
# lower : vector of lower limits for param
# identical_slopes : same slope for the dose-response curves of each substance ? (bolean)
# error_type : Residuals distribution (character chain)
#   - "Normal"
#   - "Poisson"
#   - "NB" : Negative Binomial

## 4.1. With CA as reference ----

# FONCTION CALCULATING THE SURFACE DOSE_RESPONSE KNOWING THE DOSE-RESPONSE 
# CURVES AND THE INTERACTION PARAMETERS FROM A GIVEN TYPE OF INTERACTION MODEL
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

# RSS CALCULATION OF A GIVEN DOSE-RESPONSE SURFACE 
# (NORMAL DISTRIBUTION FOR RESIDUALS)
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

# NEGATIVE LOG-LIKELIHOOD CALCULATION OF A GIVEN DOSE-RESPONSE SURFACE 
# (POISSON DISTRIBUTION FOR RESIDUALS)
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

  # Avoid 0 for numeric stability
  if (any(res_CA <= 0)) {
    return(-Inf)
  }
  
  # Negative log-likelihood
  res_CA_Poisson <- -sum(dpois(Response, lambda = res_CA, log = TRUE))
  
  return(res_CA_Poisson)
}

# NEGATIVE PSEUDO-LOG-LIKELIHOOD CALCULATION OF A GIVEN DOSE-RESPONSE SURFACE 
# (QUASI-POISSON DISTRIBUTION FOR RESIDUALS)
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
  
  # Avoid 0 for numeric stability
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
  
  # Negative Pseudo Log-likelihood
  res_CA_QuasiPoisson <- -sum(Response * log(res_CA) - res_CA - lfactorial(Response)) / phi # - Loglikelihood
  
  return(res_CA_QuasiPoisson)
}

# NEGATIVE LOG-LIKELIHOOD CALCULATION OF A GIVEN DOSE-RESPONSE SURFACE 
# (NEGATIVE BINOMIAL DISTRIBUTION FOR RESIDUALS)
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
  
  # Avoid 0 for numeric stability
  if (any(res_CA <= 0)) {
    return(-Inf)
  }
  
  neg_loglik_for_theta <- function(log_theta, y, mu) {
    theta <- exp(log_theta)  # constraint >0
    -sum(dnbinom(y, size = theta, mu = res_CA, log = TRUE))
  }
  res <- optim(par = log(1), fn = neg_loglik_for_theta, y = Response, mu = res_CA,
               method = "Brent", lower = log(1e-8), upper = log(1e6)) 
  
  est_theta <- exp(res$par)
  print(est_theta)
  
  # Negative log-likelihood
  res_CA_NB <- -sum(dnbinom(Response, size = est_theta, mu = res_CA, log = TRUE))[1]
  
  return(res_CA_NB)
}

# FIT INTERACTION PARAMETERS FOR A GIVEN TYPE OF INTERACTION MODEL
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

# FIT INTERACTION PARAMETERS FOR A GIVEN TYPE OF INTERACTION MODEL
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


# 4.2. With IA as a reference ----

# FONCTION CALCULATING THE SURFACE DOSE_RESPONSE KNOWING THE DOSE-RESPONSE 
# CURVES AND THE INTERACTION PARAMETERS FROM A GIVEN TYPE OF INTERACTION MODEL
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

# RSS CALCULATION OF A GIVEN DOSE-RESPONSE SURFACE 
# (NORMAL DISTRIBUTION FOR RESIDUALS)
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

# NEGATIVE LOG-LIKELIHOOD CALCULATION OF A GIVEN DOSE-RESPONSE SURFACE 
# (POISSON DISTRIBUTION FOR RESIDUALS)
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
  
  # Avoid 0 for numeric stability
  if (any(res_IA <= 0)) {
    return(-Inf)
  }
  
  res_IA_Poisson <- -sum(dpois(Response, lambda = res_IA, log = TRUE)) # negative Loglikelihood
  
  return(res_IA_Poisson)
}

# NEGATIVE PSEUDO-LOG-LIKELIHOOD CALCULATION OF A GIVEN DOSE-RESPONSE SURFACE 
# (QUASI-POISSON DISTRIBUTION FOR RESIDUALS)
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
  
  # Avoid 0 for numeric stability
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
  
  # Negative Pseudo Log-likelihood
  res_IA_QuasiPoisson <- -sum(Response * log(res_IA) - res_IA - lfactorial(Response)) / phi 
  
  return(res_IA_QuasiPoisson)
}

# NEGATIVE LOG-LIKELIHOOD CALCULATION OF A GIVEN DOSE-RESPONSE SURFACE 
# (NEGATIVE BINOMIAL DISTRIBUTION FOR RESIDUALS)
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
  
  # Avoid 0 for numeric stability
  if (any(res_IA <= 0)) {
    return(-Inf)
  }
  
  neg_loglik_for_theta <- function(log_theta, y, mu) {
    theta <- exp(log_theta)  # constraint >0
    -sum(dnbinom(y, size = theta, mu = res_IA, log = TRUE))
  }
  res <- optim(par = log(1), fn = neg_loglik_for_theta, y = Response, mu = res_IA,
               method = "Brent", lower = log(1e-8), upper = log(1e6)) 
  
  est_theta <- exp(res$par)
  print(est_theta)
  
  # Nagative Log-likelihood
  res_IA_NB <- -sum(dnbinom(Response, size = est_theta, mu = res_IA, log = TRUE))
  
  return(res_IA_NB)
}

# FIT INTERACTION PARAMETERS FOR A GIVEN TYPE OF INTERACTION MODEL
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

# CI CALCULATION FOR A GIVEN ESTIMATED DOSE-RESPONSE SURFACE 
# (RESAMPLING AT THE WHOLE DESIGN LEVEL, not used here)
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

# CI CALCULATION FOR A GIVEN ESTIMATED DOSE-RESPONSE SURFACE
# (RESAMPLING AT THE CONDITION LEVEL)
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



