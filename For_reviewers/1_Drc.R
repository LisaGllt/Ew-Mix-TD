library(here)
source(file = here::here("functions/fun.R"))
f_load_libraries_colors()

path_repro <- here::here("fig/Single_expo/Reproduction")
path_growth <- here::here("fig/Single_expo/Growth")

# 1. Effects on reproduction ----

## 1.1. Load data ----

data_expA_repro <- f_read_data_expA_repro()
df_expA_repro_tot <- data_expA_repro$df_expA_repro
df_expA_repro <- data_expA_repro$df_expA_repro_alive
df_expA_repro_mean <- data_expA_repro$df_expA_repro_alive_mean

Val_Ctrl <- 1e-4
Label_Ctrl <- "0 (Ctrl)"

## 1.2. Initial conditions ----

plot_D0 <- ggplot(
  data = subset(df_expA_repro_tot, t == 0),
  aes(
    x = Dose_f,
    y = w,
    color = Molecule,
    fill = Molecule,
    shape = Molecule
  )
) +
  geom_boxplot(alpha = 0.2) +
  geom_dotplot(
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.7,
    alpha = 0.7
  ) +
  facet_wrap(
    ~Molecule,
    ncol = 1,
    scales = "free"
  ) +
  scale_color_manual(values = col_molec) +
  scale_fill_manual(values = col_molec) +
  scale_shape_manual(values = shape_molec) +
  labs(x = "Dose (mg/kg)", y = "Weight (mg)") +
  theme_classic(14) +
  theme(
    legend.position = "none",
    legend.title = element_text(
      size = sizetitle - 2,
      face = "bold"
    ),
    legend.text = element_text(
      size = sizetitle - 2,
      face = "plain"
    ),
    title = element_text(
      size = sizetitle,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = sizetitle,
      face = "plain"
    ),
    axis.title.y = element_text(
      size = sizetitle,
      face = "plain"
    ),
    strip.background = element_rect(
      color = "black",
      fill = Nord_snow[3],
      size = 0,
      linetype = "solid"
    )
  )
plot_D0


# EPX
mod_aov <- aov(w ~ Dose_f, data = subset(df_expA_repro_tot, Molecule == "EPX" & t == 0))
summary(mod_aov)

dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

# IMD
mod_aov <- aov(w ~ Dose_f, data = subset(df_expA_repro_tot, Molecule == "IMD" & t == 0))
summary(mod_aov)

dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

## 1.3. Dose-response ----

drc.mod.r <- drm(
  Nb_cocoons ~ Dose,
  curveid = Molecule, # <1>
  data = df_expA_repro, # <2>
  type = "Poisson", # <3>
  fct = LL.4(
    names = c("slope", "Ymin", "Ymax", "EC50"), # <4>
    fixed = c(NA, 0, NA, NA) # <5>
  ),
  pmodels = data.frame(Molecule, Molecule, Molecule), # <6>
)

summary(drc.mod.r)

# Confidence interval

Dose_x <- expand.grid(
  exp(
    seq(
      log(Val_Ctrl), log(5000),
      by = (log(5000) - log(Val_Ctrl)) / 100
    )
  )
)

CI.r <- data.frame(
  Dose = c(Dose_x$Var1, Dose_x$Var1),
  Molecule = c(rep("EPX", length(Dose_x$Var1)), rep("IMD", length(Dose_x$Var1)))
)

pm.r <- predict(drc.mod.r, newdata = CI.r, interval = "confidence")
CI.r$p <- pm.r[, 1]
CI.r$pmin <- pm.r[, 2]
CI.r$pmax <- pm.r[, 3]

est <- coef(drc.mod.r)
ci <- confint(drc.mod.r)
params_tab <- data.frame(
  estimate  = unname(est),
  CI_lower  = ci[, 1],
  CI_upper  = ci[, 2]
)

params_tab |>
  datatable(
    options = list(
      dom = "t",
      autoWidth = TRUE,
      columnDefs = list(
        list(className = "dt-left", targets = "_all")
      )
    ),
    class = "hover"
  ) |>
  formatSignif(columns = c("estimate", "CI_lower", "CI_upper"), c(3, 3, 3))

plot_drc_r <- ggplot() +
  geom_ribbon_interactive(
    data = CI.r,
    aes(
      x = Dose,
      y = p,
      ymin = pmin, ymax = pmax,
      group = Molecule,
      fill = Molecule,
      tooltip = Molecule,
      data_id = Molecule,
    ),
    alpha = 0.2
  ) +
  geom_line_interactive(
    data = CI.r,
    aes(
      x = Dose,
      y = p,
      color = Molecule,
      tooltip = Molecule,
      data_id = Molecule,
      color = Molecule,
      shape = Molecule
    ),
    linewidth = 1.2
  ) +
  geom_point_interactive(
    data = df_expA_repro,
    aes(
      x = Dose_plot,
      y = Nb_cocoons,
      shape = Molecule,
      tooltip = Molecule,
      data_id = Molecule,
      color = Molecule,
      shape = Molecule
    ),
    cex = 2,
    alpha = 0.4
  ) +
  geom_point_interactive(
    data = df_expA_repro_mean,
    aes(
      x = Dose_plot,
      y = Nb_cocoons,
      shape = Molecule,
      tooltip = Molecule,
      data_id = Molecule,
      color = Molecule,
      shape = Molecule
    ),
    alpha = 1,
    size = 3
  ) +
  scale_color_manual(name = "", values = col_molec) +
  scale_shape_manual(name = "", values = shape_molec) +
  scale_fill_manual(name = "", values = col_molec) +
  theme_minimal() +
  scale_x_log10(
    breaks = c(Val_Ctrl, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
    labels = c(Label_Ctrl, "0.001", "0.01", "0.1", "1", "10", "100", "1000")
  ) +
  labs(
    x = "Dose (mg/kg)",
    y = "Number of cocoons produced per cosm (#)",
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(
      size = sizetitle - 2,
      face = "bold"
    ),
    title = element_text(
      size = sizetitle,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = sizetitle,
      face = "plain"
    ),
    axis.title.y = element_text(
      size = sizetitle,
      face = "plain"
    )
  )

tooltip_css <- "
  border-radius: 12px;
  color: #333;
  background-color: white;
  padding: 10px;
  font-size: 14px;
  transition: all 0.5s ease-out;
"

hover_css <- "
  filter: brightness(75%);
  cursor: pointer;
  transition: all 0.5s ease-out;
  filter: brightness(1.15);
"

interactive_plot <- girafe(ggobj = plot_drc_r) %>%
  girafe_options(opts_sizing(rescale = FALSE))

# Add interactivity
interactive_plot <- interactive_plot %>%
  girafe_options(
    opts_hover(css = hover_css),
    opts_tooltip(css = tooltip_css),
    opts_hover_inv(css = "opacity:0.3; transition: all 0.2s ease-out;")
  )

interactive_plot

plot_drc_r_f <- plot_drc_r +
  labs(title = "A. Effects on cocoon production")

## 1.4. Check model ----

# Anova model
df_expA_repro_anova <- df_expA_repro |>
  mutate(
    Condition_f = as.factor(paste0(Molecule, Dose_plot)),
    Molecule = as.factor(Molecule)
  )

anova.drc.r <- glm(
  Nb_cocoons ~ Condition_f,
  family = poisson(link = "log"),
  data = df_expA_repro_anova
)
# check_model(glm_anova)
# predict(glm_anova)

# Log-Likelihood
LL_drc <- as.numeric(logLik(drc.mod.r))
LL_anova <- as.numeric(logLik(anova.drc.r))

# Number of parameters
k_drc <- 6
k_anova <- attr(logLik(anova.drc.r), "df")

# Log-Likelihood ratio test
# H0 : the drc is OK
# H1 : the ANOVA modèle (saturated) is better
D <- -2 * (LL_drc - LL_anova)
df_test <- k_anova - k_drc

# Verification
p_val <- pchisq(D, df = df_test, lower.tail = FALSE)
cat("Lack-of-fit test: D =", round(D, 3), ", df =", df_test, ", p =", round(p_val, 4), "\n")

col_diag <- Nord_polar[4]

plot_QQ <- ggplot(mapping = aes(sample = residuals(drc.mod.r))) +
  stat_qq(alpha = 0.3) +
  stat_qq_line(color = col_diag) +
  labs(
    x = "Normal quantiles",
    y = "Residuals"
  ) +
  theme_bw()

plot_resfit <- ggplot() +
  geom_point(
    mapping = aes(
      x = fitted(drc.mod.r),
      y = residuals(drc.mod.r)
    ),
    color = col_diag,
    alpha = 0.3
  ) +
  scale_x_log10() +
  labs(
    x = "Fitted",
    y = "Residuals"
  ) +
  theme_bw()

plot <- plot_resfit + plot_QQ + plot_layout(ncol = 2)
plot

## 1.5. NOEC ----

df_expA_repro$Dose_f <- as.factor(df_expA_repro$Dose)
df_expA_repro$Dose_f <- relevel(df_expA_repro$Dose_f, ref = "0")

cat("IMD - Determination of NOEC and LOEC\n")
mod_aov <- aov(Nb_cocoons ~ Dose_f, data = subset(df_expA_repro, Molecule == "IMD"))
summary(mod_aov)
dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

cat("EPX - Determination of NOEC and LOEC\n")
mod_aov <- aov(Nb_cocoons ~ Dose_f, data = subset(df_expA_repro, Molecule == "EPX"))
summary(mod_aov)
dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

## 1.6. Effect on adult growth ----

nb_dose_EPX <- n_distinct(subset(df_expA_repro, Molecule == "EPX")$Dose)
nb_dose_IMD <- n_distinct(subset(df_expA_repro, Molecule == "IMD")$Dose)

pal_dose_EPX <- nord(
  palette = "aurora", nb_dose_EPX + 2,
  alpha = 1, reverse = T
)[-c(1, 2)]
pal_dose_IMD <- nord(
  palette = "aurora", nb_dose_IMD + 2,
  alpha = 1, reverse = T
)[-c(1, 2)]

plot_growth_EPX <- ggplot(
  data = subset(df_expA_repro, Molecule == "EPX"),
  aes(
    x = t,
    y = w,
    group = ID_cosm,
    color = as.factor(Dose)
  )
) +
  geom_point(shape = shape_EPX) +
  scale_color_manual(
    name = "Dose (mg/kg)",
    values = pal_dose_EPX
  ) +
  geom_line() +
  facet_wrap(~Dose, ncol = 3) +
  labs(
    x = "Time (d)",
    y = "Weight (mg)",
    title = "Growth of earthworms exposed to epoxiconazole"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = sizetitle - 2, face = "bold"),
    strip.background = element_rect(
      colour = "black", fill = "white",
      linetype = "solid"
    ),
    title = element_text(size = sizetitle, face = "bold"),
    axis.title.x = element_text(size = sizetitle, face = "plain"),
    axis.title.y = element_text(size = sizetitle, face = "plain")
  )

plot_growth_IMD <- ggplot(
  data = subset(df_expA_repro, Molecule == "IMD"),
  aes(
    x = t,
    y = w,
    group = ID_cosm,
    color = as.factor(Dose)
  )
) +
  geom_point(shape = shape_IMD) +
  geom_line() +
  facet_wrap(~Dose, ncol = 3) +
  scale_color_manual(
    name = "Dose (mg/kg)",
    values = pal_dose_IMD
  ) +
  labs(
    x = "Time (d)",
    y = "Weight (mg)",
    title = "Growth of earthworms exposed to imidacloprid"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = sizetitle - 2, face = "bold"),
    strip.background = element_rect(
      colour = "black", fill = "white",
      linetype = "solid"
    ),
    title = element_text(size = sizetitle, face = "bold"),
    axis.title.x = element_text(size = sizetitle, face = "plain"),
    axis.title.y = element_text(size = sizetitle, face = "plain")
  )

plot <- plot_growth_EPX + plot_growth_IMD + plot_layout(ncol = 1, heights = c(1, 4 / 3))
plot

df_expA_repro_grate <- df_expA_repro |>
  dplyr::select(ID_cosm, t, L, Dose, Dose_plot, Molecule) |>
  pivot_wider(names_from = t, values_from = L, names_prefix = "L_t") |>
  mutate(a_i = (L_t28 - L_t0) / 28)

df_expA_repro_grate_mean <- df_expA_repro_grate |>
  aggregate(a_i ~ Molecule + Dose + Dose_plot, FUN = "mean")

p <- ggplot() +
  geom_point(
    data = df_expA_repro_grate,
    aes(
      x = Dose_plot,
      y = a_i,
      color = Molecule,
      shape = Molecule
    ),
    alpha = 0.5,
  ) +
  geom_point(
    data = df_expA_repro_grate_mean,
    aes(
      x = Dose_plot,
      y = a_i,
      color = Molecule,
      shape = Molecule
    ),
    alpha = 1,
    size = 3
  ) +
  scale_x_log10(
    breaks = c(Val_Ctrl, 0.01, 0.1, 1, 10, 100), # inclure la valeur fictive
    labels = c(Label_Ctrl, "0.01", "0.1", "1", "10", "100")
  ) +
  facet_wrap(~Molecule, scales = "free_x") +
  scale_color_manual(values = col_molec) +
  scale_fill_manual(values = col_molec) +
  scale_shape_manual(values = shape_molec) +
  labs(
    x = "Dose (mg/kg)",
    y = "Growth rate" # ,
    # title = "Effect of the pesticides on adult growth rates"
  ) +
  ylim(NA, 0.05) +
  theme_minimal(base_size = 15) + # Thème minimaliste
  theme(
    plot.title = element_text(
      # hjust = 0.5,
      face = "bold"
    ),
    legend.position = "none",
    strip.background = element_rect(
      color = "black",
      fill = Nord_snow[3],
      size = 0,
      linetype = "solid"
    ),
    strip.text = element_text(face = "bold")
  )
p

# NOEC

df_expA_repro_grate$Dose_f <- as.factor(df_expA_repro_grate$Dose)
df_expA_repro_grate$Dose_f <- relevel(df_expA_repro_grate$Dose_f, ref = "0")

cat("IMD - Determination of NOEC and LOEC\n")
mod_aov <- aov(a_i ~ Dose_f, data = subset(df_expA_repro_grate, Molecule == "IMD"))
summary(mod_aov)
dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

cat("EPX - Determination of NOEC and LOEC\n")
mod_aov <- aov(a_i ~ Dose_f, data = subset(df_expA_repro_grate, Molecule == "EPX"))
summary(mod_aov)
dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

# 2. Effect on juvenile growth ----

## 2.1. Load data ----

data_expA_growth <- f_read_data_expA_growth()
df_expA_growth_tot <- data_expA_growth$df_expA_growth_tot
df_expA_growth <- data_expA_growth$df_expA_growth

## 2.2. Initial state ----

plot_D0 <- ggplot(
  data = subset(df_expA_growth_tot, t == 0),
  aes(
    x = Dose_f,
    y = w,
    color = Molecule,
    fill = Molecule
  )
) +
  geom_boxplot(alpha = 0.2) +
  geom_dotplot(
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.7,
    alpha = 0.7
  ) +
  labs(
    x = "Dose (mg/kg)",
    y = "Weight (mg)"
  ) +
  facet_wrap(
    ~Molecule,
    ncol = 1,
    scales = "free"
  ) +
  scale_color_manual(values = col_molec) +
  scale_fill_manual(values = col_molec) +
  theme_classic(14) +
  theme(
    legend.position = "none",
    legend.title = element_text(
      size = sizetitle - 2,
      face = "bold"
    ),
    legend.text = element_text(
      size = sizetitle - 2,
      face = "plain"
    ),
    title = element_text(
      size = sizetitle,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = sizetitle,
      face = "plain"
    ),
    axis.title.y = element_text(
      size = sizetitle,
      face = "plain"
    ),
    strip.background = element_rect(
      color = "black",
      fill = "grey95",
      size = 0,
      linetype = "solid"
    )
  )
plot_D0

# ANOVA
mod_aov <- aov(w ~ Dose_f, data = subset(df_expA_growth_tot, Molecule == "EPX" & t == 0))
summary(mod_aov)

dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

# ANOVA
mod_aov <- aov(w ~ Dose_f, data = subset(df_expA_growth_tot, Molecule == "IMD" & t == 0))
summary(mod_aov)

dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

## 2.3. Growth curves ----

nb_dose_EPX <- n_distinct(subset(df_expA_growth_tot, Molecule == "EPX")$Dose)
nb_dose_IMD <- n_distinct(subset(df_expA_growth_tot, Molecule == "IMD")$Dose)

pal_dose_EPX <- nord(
  palette = "aurora",
  nb_dose_EPX + 2,
  alpha = 1,
  reverse = T
)[-c(1, 2)]
pal_dose_IMD <- nord(
  palette = "aurora",
  nb_dose_IMD + 2,
  alpha = 1,
  reverse = T
)[-c(1, 2)]

plot_growth_EPX <- ggplot(
  data = subset(df_expA_growth_tot, Molecule == "EPX"),
  aes(
    x = t,
    y = w,
    group = ID,
    color = as.factor(Dose)
  ),
  shape = shape_EPX
) +
  geom_point() +
  scale_color_manual(
    name = "Dose (mg/kg)",
    values = pal_dose_EPX
  ) +
  geom_line() +
  facet_wrap(
    ~Dose,
    ncol = 3
  ) +
  labs(
    x = "Time (d)",
    y = "Weight (mg)",
    title = "Epoxiconazole"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(
      size = sizetitle - 2,
      face = "bold"
    ),
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linetype = "solid"
    ),
    title = element_text(
      size = sizetitle,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = sizetitle,
      face = "plain"
    ),
    axis.title.y = element_text(
      size = sizetitle,
      face = "plain"
    )
  )

plot_growth_IMD <- ggplot(
  data = subset(df_expA_growth_tot, Molecule == "IMD"),
  aes(
    x = t,
    y = w,
    group = ID,
    color = as.factor(Dose)
  ),
  shape = shape_IMD
) +
  geom_point() +
  geom_line() +
  facet_wrap(
    ~Dose,
    ncol = 3
  ) +
  scale_color_manual(name = "Dose (mg/kg)", values = pal_dose_IMD) +
  labs(
    x = "Time (d)",
    y = "Weight (mg)",
    title = "Imidacloprid"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(
      size = sizetitle - 2,
      face = "bold"
    ),
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linetype = "solid"
    ),
    title = element_text(
      size = sizetitle,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = sizetitle,
      face = "plain"
    ),
    axis.title.y = element_text(
      size = sizetitle,
      face = "plain"
    )
  )

plot <- plot_growth_EPX + plot_growth_IMD + plot_layout(ncol = 1, heights = c(1, 4 / 3))
plot

## 2.4. Growth rate calculations ----

f_growth_rate <- function(data) {
  a_init <- 0.1
  a_i <- c()
  
  # Calculation of the growth rate for each individual
  for (i in subset(data, t == 0)$ID) {
    data_i <- data.frame(
      L = subset(data, ID == i)$L,
      t = subset(data, ID == i)$t
    )
    L0_i <- subset(data, ID == i & t == 0)$L
    
    if (length(na.omit(data_i$L)) == 2) { # If only 2 data points
      a.Est <- (subset(data_i, !t == 0)$L - subset(data_i, t == 0)$L) /
        subset(data_i, !t == 0)$t
    } else { # If 3 data points (or more)
      EPX.nlsHill <- nls(
        L ~ L0_i + t * a,
        start = list(a = a_init),
        trace = T,
        data = data_i
      )
      a.Est <- EPX.nlsHill$m$getPars()[1]
    }
    a_i <- rbind(a_i, a.Est)
  }
  
  return(a_i)
}

# Growth rate calculation for individuals alive at the end of the experiment
a_i <- f_growth_rate(df_expA_growth) |> as.data.frame()

df_rates_alive <- cbind(subset(df_expA_growth, t == 0), a_i) |>
  as.data.frame()

df_rates_controls <- subset(df_rates_alive, Molec %in% c("Ctrl_1", "Ctrl_2")) |>
  aggregate(a ~ Molec, FUN = mean)

df_rates_alive <- df_rates_alive |>
  mutate(
    r_a = case_when(
      Lot == "A" ~ a / subset(df_rates_controls, Molec == "Ctrl_1")$a,
      Lot == "B" ~ a / subset(df_rates_controls, Molec == "Ctrl_2")$a
    ) * 100
  )

df_rates_alive_mean <- df_rates_alive |>
  group_by(Dose, Molecule, Dose_plot) |>
  summarise(
    r_a = mean(r_a, na.rm = TRUE),
    a = mean(a, na.rm = TRUE),
    .groups = "drop"
  )

## 2.5 Dose-response curves ----

drc.mod.g <- drm(
  a ~ Dose, Molecule,
  data = df_rates_alive,
  fct = LL.4(
    names = c("slope", "Ymin", "Ymax", "EC50"),
    fixed = c(NA, NA, NA, NA)
  ),
  pmodels = data.frame(Molecule, 1, 1, Molecule),
  lowerl = c(-Inf, -Inf, -Inf, 1e-4)
)
summary(drc.mod.g)

Dose_x <- expand.grid(
  exp(
    seq(
      log(0.000001), log(5000),
      by = (log(5000) - log(0.000001)) / 100
    )
  )
)

CI.g <- data.frame(
  Dose = c(Dose_x$Var1, Dose_x$Var1),
  Molecule = c(rep("EPX", length(Dose_x$Var1)), rep("IMD", length(Dose_x$Var1)))
)

pm.g <- predict(drc.mod.g, newdata = CI.g, interval = "confidence")
CI.g$p <- pm.g[, 1]
CI.g$pmin <- pm.g[, 2]
CI.g$pmax <- pm.g[, 3]

plot_drc_g <- ggplot() +
  # geom_ribbon_interactive(
  #   data=CI.g,
  #   aes(
  #     x=Dose,
  #     y=p,
  #     ymin=pmin, ymax=pmax,
  #     group=Molecule,
  #     fill=Molecule,
  #     tooltip = Molecule,
  #     data_id = Molecule,
  #     ),
  #   alpha=0.2)+
  # geom_line_interactive(
  #   data = CI.g,
  #   aes(
  #     x=Dose,
  #     y=p,
  #     color=Molecule,
  #     tooltip = Molecule,
  #     data_id = Molecule,
  #     color=Molecule,
  #     shape = Molecule
  #     ),
  #   linewidth = 1.2
  #   )+
  geom_point_interactive(
    data = df_rates_alive,
    aes(
      x = Dose_plot,
      y = a,
      shape = Molecule,
      tooltip = Molecule,
      data_id = Molecule,
      color = Molecule,
      shape = Molecule
    ),
    cex = 2,
    alpha = 0.4
  ) +
  geom_point_interactive(
    data = df_rates_alive_mean,
    aes(
      x = Dose_plot,
      y = a,
      shape = Molecule,
      tooltip = Molecule,
      data_id = Molecule,
      color = Molecule,
      shape = Molecule
    ),
    alpha = 1,
    size = 3
  ) +
  scale_color_manual(name = "", values = col_molec) +
  scale_shape_manual(name = "", values = shape_molec) +
  scale_fill_manual(name = "", values = col_molec) +
  theme_minimal() +
  scale_x_log10(
    breaks = c(Val_Ctrl, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
    labels = c(Label_Ctrl, "0.001", "0.01", "0.1", "1", "10", "100", "1000")
  ) +
  labs(
    x = "Dose (mg/kg)",
    y = expression(paste("Growth rate (", mg^{
      1 / 3
    }, "/day)")),
    # y = "Relative growth rate (%)",
    # title = "B. Effects on juvenile growth"
  ) +
  theme(
    legend.position = "none",
    legend.title = element_text(
      size = sizetitle - 2,
      face = "bold"
    ),
    title = element_text(
      size = sizetitle,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = sizetitle,
      face = "plain"
    ),
    axis.title.y = element_text(
      size = sizetitle,
      face = "plain"
    )
  )

tooltip_css <- "
  border-radius: 12px;
  color: #333;
  background-color: white;
  padding: 10px;
  font-size: 14px;
  transition: all 0.5s ease-out;
"

hover_css <- "
  filter: brightness(75%);
  cursor: pointer;
  transition: all 0.5s ease-out;
  filter: brightness(1.15);
"

interactive_plot <- girafe(ggobj = plot_drc_g) %>%
  girafe_options(opts_sizing(rescale = FALSE))

# Add interactivity
interactive_plot <- interactive_plot %>%
  girafe_options(
    opts_hover(css = hover_css),
    opts_tooltip(css = tooltip_css),
    opts_hover_inv(css = "opacity:0.3; transition: all 0.2s ease-out;")
  )

plot_drc_g_f <- plot_drc_g +
  labs(title = "B. Effects on juvenile growth")

## 2.6. Check model ----

# Anova model
df_rates_alive <- df_rates_alive |>
  mutate(
    Condition_f = as.factor(paste0(Molecule, Dose_plot)),
    Molecule = as.factor(Molecule)
  )

anova.drc.g <- glm(
  a ~ Condition_f,
  data = df_rates_alive
)
# check_model(glm_anova)
# predict(glm_anova)

# Log-likelihood
LL_drc <- as.numeric(logLik(drc.mod.g))
LL_anova <- as.numeric(logLik(anova.drc.g))

# Number of parameters
k_drc <- 6
k_anova <- attr(logLik(anova.drc.g), "df")

# Log-likelihood ratio test
# H0 : the drc model is OK
# H1 : the ANOVA model (saturated) is better
D <- -2 * (LL_drc - LL_anova)
df_test <- k_anova - k_drc

# Vérification

p_val <- pchisq(D, df = df_test, lower.tail = FALSE)
cat("Lack-of-fit test: D =", round(D, 3), ", df =", df_test, ", p =", round(p_val, 4), "\n")

col_diag <- Nord_polar[4]

plot_QQ <- ggplot(
  mapping = aes(
    sample = residuals(drc.mod.g)
  )
) +
  stat_qq(alpha = 0.3) +
  stat_qq_line(color = col_diag) +
  labs(
    x = "Normal quantiles",
    y = "Residuals"
  ) +
  theme_minimal()

plot_resfit <- ggplot() +
  geom_point(
    mapping = aes(
      x = fitted(drc.mod.g),
      y = residuals(drc.mod.g)
    ),
    color = col_diag,
    alpha = 0.3
  ) +
  labs(
    x = "Fitted",
    y = "Residuals"
  ) +
  scale_x_log10() +
  theme_minimal()

plot <- plot_resfit + plot_QQ + plot_layout(ncol = 2)
plot

## 2.7. NOEC ----

# ANOVA
cat("IMD - Determination of NOEC and LOEC\n")
mod_aov <- aov(a ~ Dose_f, data = subset(df_rates_alive, Molecule == "IMD"))
summary(mod_aov)

dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

cat("EPX - Determination of NOEC and LOEC\n")
mod_aov <- aov(a ~ Dose_f, data = subset(df_rates_alive, Molecule == "EPX"))
summary(mod_aov)

dunnett_test <- glht(mod_aov, linfct = mcp(Dose_f = "Dunnett"))
summary(dunnett_test)

## 2.8. Survival ----

df_dead_g <- df_expA_growth_tot |>
  filter(t == 28) |>
  group_by(Dose_f, Molecule) |>
  summarise(
    n_D = sum(Status == "D", na.rm = TRUE), # compte lignes avec Status "D"
    n_total = n(), # total lignes par condition
    .groups = "drop"
  )

p <- ggplot(
  data = df_dead_g,
  aes(
    x = Dose_f,
    y = n_D,
    color = Molecule,
    fill = Molecule
  )
) +
  geom_col(
    alpha = 0.5
  ) +
  geom_text(
    aes(label = round(n_D)),
    vjust = -0.5,
    size = 4
  ) +
  labs(
    x = "Dose (mg/kg)",
    y = "Number of dead individuals (#)"
  ) +
  facet_wrap(
    ~Molecule,
    ncol = 1,
    scales = "free"
  ) +
  scale_color_manual(values = col_molec) +
  scale_fill_manual(values = col_molec) +
  lims(y = c(0, 4.5)) +
  theme_classic(12) +
  theme(
    legend.position = "none",
    legend.title = element_text(
      size = sizetitle - 2,
      face = "bold"
    ),
    legend.text = element_text(
      size = sizetitle - 2,
      face = "plain"
    ),
    title = element_text(
      size = sizetitle,
      face = "bold"
    ),
    strip.background = element_rect(
      color = "black",
      fill = "grey95",
      size = 0,
      linetype = "solid"
    )
  )
p


