library(tidyverse)
library(readxl)
source(file="fun.R")

N_draws <- 5000

load(here::here("mod/Docs_for_CI_calc/Data_mix_CI.RData"))

CI_CASA <- f_Jonker_CI_calculations_condition(
  df_data = df_expB_repro_Jonker,
  N_draws = N_draws,
  interact = "SA",
  reference = "CA",
  error_type = "Poisson"
)
save(CI_CASA, file = "CI_cond_CASA.RData")

CI_CADL <- f_Jonker_CI_calculations_condition(
  df_data = df_expB_repro_Jonker,
  N_draws = N_draws,
  interact = "DL",
  reference = "CA",
  error_type = "Poisson"
)

save(CI_CADL, file = "CI_cond_CADL.RData")

CI_CADR <- f_Jonker_CI_calculations_condition(
  df_data = df_expB_repro_Jonker,
  N_draws = N_draws,
  interact = "DR",
  reference = "CA",
  error_type = "Poisson"
)

save(CI_CADR, file = "CI_cond_CADR.RData")


CI_IASA <- f_Jonker_CI_calculations_condition(
  df_data = df_expB_repro_Jonker,
  N_draws = N_draws,
  interact = "SA",
  reference = "IA",
  error_type = "Poisson"
)
save(CI_IASA, file = "CI_cond_IASA.RData")

CI_IADL <- f_Jonker_CI_calculations_condition(
  df_data = df_expB_repro_Jonker,
  N_draws = N_draws,
  interact = "DL",
  reference = "IA",
  error_type = "Poisson"
)

save(CI_IADL, file = "CI_cond_IADL.RData")

CI_IADR <- f_Jonker_CI_calculations_condition(
  df_data = df_expB_repro_Jonker,
  N_draws = N_draws,
  interact = "DR",
  reference = "IA",
  error_type = "Poisson"
)

save(CI_IADR, file = "CI_cond_IADR.RData")