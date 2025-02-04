Simulation: Results
================
<br> Updated on April 17, 2024

**Overview**: This code looks into simulation study results and creates
visualizations for paper.

Install correct version of `listenr`:

``` r
# install.packages(
#   "../listenr-v0.3.4.tar.gz",
#   repos = NULL,
#   type = "source"
# )
```

Load packages:

``` r
library(dplyr)
library(ggplot2)
library(listenr)
library(tidyr)
library(viridis)
```

# Example Simulated Datasets

Specify values simulating datasets:

``` r
# Number of time points and spatial locations
nT = 70
nx = 10
ny = 10

# Parameter values for data simulation
x1mean = 20
x2mean = 45
xsd = 6

# Sigma parameter values
grid <-
  data.frame(
    sigmaX = c(0.2,4),
    sigmaDelta = c(0.2,4),
    sigmaEpsilon = c(0.2,4)
  )
```

Generate data:

``` r
out_list <- list()
for (i in 1:nrow(grid)) {
  dat = simulate_synthetic_data(
    nT = nT,
    nx = nx,
    ny = ny,
    x1mean = x1mean,
    x2mean = x2mean,
    xsd = xsd,
    sigmaX = grid$sigmaX[i],
    sigmaDelta = grid$sigmaDelta[i],
    sigmaEpsilon = grid$sigmaEpsilon[i],
    phiX = .5,
    phiDelta = .5,
    rhoX = .5,
    rhoDelta = .5,
    seed = 12345
  )
  dat$sigmaX <- grid$sigmaX[i]
  dat$sigmaDelta <- grid$sigmaDelta[i]
  dat$sigmaEpsilon <- grid$sigmaEpsilon[i]
  dat$group <- i
  out_list[[i]] <- dat
}
outdf <- do.call("rbind", out_list)
```

Clean up generated data:

``` r
moutdf <-
  outdf |>
  filter(time <= 70) |>
  mutate(temp = factor(sigmaX, labels = c(expression(
    paste(sigma[z], "=", sigma[delta], "=", sigma[epsilon], "= 0.2")
  ), expression(
    paste(sigma[z], "=", sigma[delta], "=", sigma[epsilon], "= 4")
  )))) |>
  group_by(time, temp, group) |>
  summarize(mval = mean(val),
            mx1 = mean(x1),
            mx2 = mean(x2)) |>
  pivot_longer(mval:mx2, names_to = "variable", values_to = "value")
```

Create color palette:

``` r
cbPalette <-
  c(
    "#999999",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )
```

Visualize examples of spatially averaged data:

``` r
# Create/load plot
fp_p_sim_means = "../figs/simulation_spatially_averaged.png"
if (!file.exists(fp_p_sim_means)) {
  p_sim_means <-
    moutdf |>
    ggplot(aes(
      x = time,
      y = value,
      group = interaction(group, variable),
      color = variable
    )) +
    geom_line(linewidth = 0.75) +
    facet_grid(. ~ temp, labeller = label_parsed) +
    scale_colour_manual(
      values = c("black", cbPalette[2:3]),
      labels = c(expression(Z[y]), expression(Z[1]), expression(Z[2]))
    ) +
    theme_bw(base_size = 24) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 24),
      legend.title = element_text(size = 24),
      legend.text = element_text(size = 24),
      axis.title = element_text(size = 24)
    ) +
    labs(
      x = "Time",
      y = "Value",
      color = "Variable",
      title = expression(
        paste("Spatially averaged values of variables ", 
              Z[y], ", ", Z[1], ", ", Z[2])
      )
    )
  ggsave(
    filename = fp_p_sim_means,
    plot = p_sim_means,
    dpi = 500,
    height = 6.5,
    width = 13
  )
}
knitr::include_graphics(fp_p_sim_means)
```

<img src="../figs/simulation_spatially_averaged.png" width="6500" />

Visualize spatial example of data:

``` r
# Create/load plot
fp_p_sim_spat = "../figs/simulation_heatmaps.png"
if (!file.exists(fp_p_sim_spat)) {
  p_sim_spat <-
    dat |>
    filter(time <= 70) |>
    ggplot(aes(x = easting, y = northing, fill = val)) +
    geom_tile() +
    facet_wrap(~ time, ncol = 10) +
    scale_fill_viridis(option = "magma") +
    labs(
      fill = latex2exp::TeX("$Z_{Y,t}(\\textbf{s}_i)$"),
      x = "",
      y = ""#,
      #title = "Grid ordered by times 1 to 70: left to right, top to bottom"
    ) +
    theme_bw(base_size = 26) +
    theme(
      aspect.ratio = 1,
      strip.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 20),
      axis.title = element_blank(),
      strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")),
      legend.position = "bottom",
      legend.key.width = unit(2, 'cm')
    )
  ggsave(
    filename = fp_p_sim_spat,
    plot = p_sim_spat,
    dpi = 500,
    height = 12,
    width = 12
  )
}
knitr::include_graphics(fp_p_sim_spat)
```

<img src="../figs/simulation_heatmaps.png" width="6000" />

# Simulation Results

Read in results:

``` r
res_sigma <- read.csv("../output/simulation_res_sigma.csv")
res_phi_rho <- read.csv("../output/simulation_res_phi_rho.csv")
```

Aggregate results:

``` r
agg_res_sigma <-
  res_sigma |>
  group_by(
    feat_imp,
    retraining,
    spatial,
    variable,
    vars_adj,
    t_adj,
    t_forecasted,
    simGroup,
    sigmaX,
    sigmaDelta,
    sigmaEpsilon,
    phiX,
    phiDelta,
    rhoX,
    rhoDelta,
    nblocks
  ) |>
  summarize(m_importance = mean(importance),
            s_importance = sd(importance)) |>
  mutate(
    phiX = paste0("phiX=", phiX),
    phiDelta = paste0("phiDelta=", phiDelta),
    rhoX = paste0("rhoX=", rhoX),
    rhoDelta = paste0("rhoDelta=", rhoDelta),
    sigmaX = paste0("sigmaX=", sigmaX),
    sigmaDelta = paste0("sigmaDelta=", sigmaDelta),
    sigmaEpsilon = paste0("sigmaEpsilon=", sigmaEpsilon)
  )

agg_res_phi_rho <- res_phi_rho |>
  group_by(
    feat_imp,
    retraining,
    spatial,
    variable,
    vars_adj,
    t_adj,
    t_forecasted,
    simGroup,
    sigmaX,
    sigmaDelta,
    sigmaEpsilon,
    phiX,
    phiDelta,
    rhoX,
    rhoDelta,
    nblocks
  ) |>
  summarize(m_importance = mean(importance),
            s_importance = sd(importance)) |>
  mutate(
    phiX = paste0("phiX=", phiX),
    phiDelta = paste0("phiDelta=", phiDelta),
    rhoX = paste0("rhoX=", rhoX),
    rhoDelta = paste0("rhoDelta=", rhoDelta),
    sigmaX = paste0("sigmaX=", sigmaX),
    sigmaDelta = paste0("sigmaDelta=", sigmaDelta),
    sigmaEpsilon = paste0("sigmaEpsilon=", sigmaEpsilon)
  )
```

## Main Paper

Plot comparing PFI and ZFI:

``` r
# Create/load plot
fp_p_sim_pfi_zfi = "../figs/simulation_zfi_pfi_comparison.png"
if (!file.exists(fp_p_sim_pfi_zfi)) {
  p_sim_pfi_zfi <-
    agg_res_sigma |>
    mutate(feat_imp2 = ifelse(feat_imp == "PFI", "Permuted", "Zeroed")) |>
    filter(
      t_adj <= 70,
      sigmaX != "sigmaX=2",
      sigmaDelta != "sigmaDelta=2",
      sigmaEpsilon == "sigmaEpsilon=4"
    ) |>
    mutate(
      sigmaX = factor(
        sigmaX,
        levels = c("sigmaX=0.2", "sigmaX=2", "sigmaX=4"),
        labels = c(expression(paste(sigma[z], "=0.2")),
                   expression(paste(sigma[z], "=2")),
                   expression(paste(sigma[z], "=4")))
      ),
      sigmaDelta = factor(
        sigmaDelta,
        levels = c("sigmaDelta=0.2", "sigmaDelta=2", "sigmaDelta=4"),
        labels = c(expression(paste(sigma[delta], "=0.2")),
                   expression(paste(sigma[delta], "=2")),
                   expression(paste(sigma[delta], "=4")))
      ),
      variable = factor(
        variable,
        levels = c("x1", "x2"),
        labels = c(expression(Z[1]), expression(Z[2]))
      )
    ) |>
    ggplot(aes(
      x = t_adj,
      y = m_importance,
      group = interaction(feat_imp, as.factor(simGroup)),
      color = (feat_imp2)
    )) +
    geom_line(alpha = 0.75) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               alpha = 0.25) +
    facet_grid(nblocks + sigmaDelta ~ variable + sigmaX ,
               scales = "free_y",
               labeller = label_parsed) +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1[2], "black")) +
    labs(
      x = "Time",
      y = "Feature importance",
      color = "",
      subtitle = ("Note: y-axis scales differ by row"),
      title = expression(paste("Feature Importances for ", sigma[epsilon], "=4"))
    ) +
    theme_bw(base_size = 18)
  ggsave(
    filename = fp_p_sim_pfi_zfi,
    plot = p_sim_pfi_zfi,
    dpi = 500,
    height = 11,
    width = 16
  )
}
knitr::include_graphics(fp_p_sim_pfi_zfi)
```

<img src="../figs/simulation_zfi_pfi_comparison.png" width="8000" />

Plot comparing number of blocks with ZFI:

``` r
# Create/load plot
fp_p_sim_blocks = "../figs/simulation_zfi_nblock.png"
if (!file.exists(fp_p_sim_blocks)) {
  p_sim_blocks <-
    agg_res_sigma |>
    mutate(feat_imp2 = ifelse(feat_imp == "PFI", "Permuted", "Zereod")) |>
    filter(
      feat_imp == "ZFI",
      sigmaX != "sigmaX=2",
      sigmaDelta != "sigmaDelta=2",
      sigmaEpsilon != "sigmaEpsilon=2"
    ) |>
    filter(
      t_adj <= 70,
      sigmaX != "sigmaX=2",
      sigmaDelta != "sigmaDelta=2",
      sigmaEpsilon != "sigmaEpsilon=1"
    ) |>
    mutate(
      sigmaX = factor(
        sigmaX,
        levels = c("sigmaX=0.2", "sigmaX=2", "sigmaX=4"),
        labels = c(expression(paste(sigma[z], "=0.2")),
                   expression(paste(sigma[z], "=2")),
                   expression(paste(sigma[z], "=4")))
      ),
      sigmaDelta = factor(
        sigmaDelta,
        levels = c("sigmaDelta=0.2", "sigmaDelta=2", "sigmaDelta=4"),
        labels = c(expression(paste(sigma[delta], "=0.2")),
                   expression(paste(sigma[delta], "=2")),
                   expression(paste(sigma[delta], "=4")))
      ),
      sigmaEpsilon = factor(
        sigmaEpsilon,
        levels = c("sigmaEpsilon=0.2", "sigmaEpsilon=2", "sigmaEpsilon=4"),
        labels = c(expression(paste(sigma[epsilon], "=0.2")),
                   expression(paste(sigma[epsilon], "=2")),
                   expression(paste(sigma[epsilon], "=4")))
      ),
      variable = factor(
        variable,
        levels = c("x1", "x2"),
        labels = c(expression(Z[1]), expression(Z[2]))
      )
    ) |>
    ggplot(aes(
      x = t_adj,
      y = m_importance,
      group = interaction(feat_imp, as.factor(simGroup)),
      color = as.factor(nblocks)
    )) +
    geom_line(alpha = 0.75) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               alpha = 0.25) +
    facet_grid(sigmaDelta + sigmaEpsilon ~ variable + sigmaX, labeller = label_parsed) +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1[3:2], "black")) +
    labs(
      x = "Time",
      y = "Feature importance",
      color = "Block size",
      title = "Zeroed Feature Importances Only"
    ) +
    theme_bw(base_size = 17)
  ggsave(
    filename = fp_p_sim_blocks,
    plot = p_sim_blocks,
    dpi = 500,
    height = 8,
    width = 16
  )
}
knitr::include_graphics(fp_p_sim_blocks)
```

<img src="../figs/simulation_zfi_nblock.png" width="8000" />

## Supplemental

Full PFI results:

``` r
# Create/load plot
fp_p_sim_full_pfi = "../figs/simulation_pfi_supplement.png"
if (!file.exists(fp_p_sim_full_pfi)) {
  p_sim_full_pfi <-
    agg_res_sigma |>
    filter(feat_imp == "PFI") |>
    filter(t_adj <= 70) |>
    mutate(
      sigmaX = factor(
        sigmaX,
        levels = c("sigmaX=0.2", "sigmaX=2", "sigmaX=4"),
        labels = c(expression(paste(sigma[z], "=0.2")),
                   expression(paste(sigma[z], "=2")),
                   expression(paste(sigma[z], "=4")))
      ),
      sigmaDelta = factor(
        sigmaDelta,
        levels = c("sigmaDelta=0.2", "sigmaDelta=2", "sigmaDelta=4"),
        labels = c(expression(paste(sigma[delta], "=0.2")),
                   expression(paste(sigma[delta], "=2")),
                   expression(paste(sigma[delta], "=4")))
      ),
      sigmaEpsilon = factor(
        sigmaEpsilon,
        levels = c("sigmaEpsilon=0.2", "sigmaEpsilon=2", "sigmaEpsilon=4"),
        labels = c(expression(paste(sigma[epsilon], "=0.2")),
                   expression(paste(sigma[epsilon], "=2")),
                   expression(paste(sigma[epsilon], "=4")))
      ),
      variable = factor(
        variable,
        levels = c("x1", "x2"),
        labels = c(expression(Z[1]), expression(Z[2]))
      )
    ) |>
    ggplot(aes(
      x = t_adj,
      y = m_importance,
      group = interaction(feat_imp, as.factor(simGroup)),
      color = as.factor(nblocks)
    )) +
    geom_line(alpha = 0.75) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               alpha = 0.25) +
    facet_grid(sigmaDelta + sigmaEpsilon ~ variable + sigmaX, labeller = label_parsed) +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1[3:2], "black")) +
    labs(
      x = "Time",
      y = "Feature importance",
      color = "Block size",
      title = "Permuted Feature Importances"
    ) +
    theme_bw(base_size = 16)
  ggsave(
    filename = fp_p_sim_full_pfi,
    plot = p_sim_full_pfi,
    dpi = 500,
    height = 14,
    width = 15
  )
}
knitr::include_graphics(fp_p_sim_full_pfi)
```

<img src="../figs/simulation_pfi_supplement.png" width="7500" />

Full ZFI results:

``` r
# Create/load plot
fp_p_sim_full_zfi = "../figs/simulation_zfi_supplement.png"
if (!file.exists(fp_p_sim_full_zfi)) {
  p_sim_full_zfi <-
    agg_res_sigma |>
    filter(feat_imp == "ZFI") |>
    filter(t_adj <= 70) |>
    mutate(
      sigmaX = factor(
        sigmaX,
        levels = c("sigmaX=0.2", "sigmaX=2", "sigmaX=4"),
        labels = c(expression(paste(sigma[z], "=0.2")),
                   expression(paste(sigma[z], "=2")),
                   expression(paste(sigma[z], "=4")))
      ),
      sigmaDelta = factor(
        sigmaDelta,
        levels = c("sigmaDelta=0.2", "sigmaDelta=2", "sigmaDelta=4"),
        labels = c(expression(paste(sigma[delta], "=0.2")),
                   expression(paste(sigma[delta], "=2")),
                   expression(paste(sigma[delta], "=4")))
      ),
      sigmaEpsilon = factor(
        sigmaEpsilon,
        levels = c("sigmaEpsilon=0.2", "sigmaEpsilon=2", "sigmaEpsilon=4"),
        labels = c(expression(paste(sigma[epsilon], "=0.2")),
                   expression(paste(sigma[epsilon], "=2")),
                   expression(paste(sigma[epsilon], "=4")))
      ),
      variable = factor(
        variable,
        levels = c("x1", "x2"),
        labels = c(expression(Z[1]), expression(Z[2]))
      )
    ) |>
    ggplot(aes(
      x = t_adj,
      y = m_importance,
      group = interaction(feat_imp, as.factor(simGroup)),
      color = as.factor(nblocks)
    )) +
    geom_line(alpha = 0.75) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               alpha = 0.25) +
    facet_grid(sigmaDelta + sigmaEpsilon ~ variable + sigmaX, labeller = label_parsed) +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1[3:2], "black")) +
    labs(
      x = "Time",
      y = "Feature importance",
      color = "Block size",
      title = "Zeroed Feature Importances",
    ) +
    theme_bw(base_size = 16)
  ggsave(
    filename = fp_p_sim_full_zfi,
    plot = p_sim_full_zfi,
    dpi = 500,
    height = 14,
    width = 15
  )
}
knitr::include_graphics(fp_p_sim_full_zfi)
```

<img src="../figs/simulation_zfi_supplement.png" width="7500" />

Comparison of PFI with varied rho and phi:

``` r
# Create/load plot
fp_p_sim_rho_pfi = "../figs/simulation_pfi_rho_phi.png"
if (!file.exists(fp_p_sim_rho_pfi)) {
  p_sim_rho_pfi <-
    agg_res_phi_rho |>
    filter(feat_imp == "PFI", nblocks == 3) |>
    filter(t_adj <= 70) |>
    mutate(
      rhoX = factor(
        rhoX,
        levels = c("rhoX=0.2", "rhoX=0.8"),
        labels = c(expression(paste(rho[Z], "=0.2")), expression(paste(rho[Z], "=0.8")))
      ),
      rhoDelta = factor(
        rhoDelta,
        levels = c("rhoDelta=0.2", "rhoDelta=0.8"),
        labels = c(expression(paste(rho[delta], "=0.2")), expression(paste(rho[delta], "=0.8")))
      ),
      phiX = factor(
        phiX,
        levels = c("phiX=0.2", "phiX=0.8"),
        labels = c(expression(paste(phi[Z], "=0.2")), expression(paste(phi[Z], "=0.8")))
      ),
      phiDelta = factor(
        phiDelta,
        levels = c("phiDelta=0.2", "phiDelta=0.8"),
        labels = c(expression(paste(phi[delta], "=0.2")), expression(paste(phi[delta], "=0.8")))
      ),
      variable = factor(
        variable,
        levels = c("x1", "x2"),
        labels = c("Z1", "Z2")
      )
    ) |>
    ggplot(aes(
      x = t_adj,
      y = m_importance,
      group = interaction(variable, as.factor(simGroup)),
      color = variable
    )) +
    geom_line(alpha = 0.75) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               alpha = 0.25) +
    facet_grid(phiX + phiDelta ~ rhoX + rhoDelta, labeller = label_parsed) +
    scale_color_manual(values = cbPalette[c(7, 6)]) +
    labs(
      x = "Time",
      y = "Feature importance",
      color = "Variable",
      title = "Permuted Feature Importances (nblocks=3)"
    ) +
    theme_bw(base_size = 18)
  ggsave(
    filename = fp_p_sim_rho_pfi,
    plot = p_sim_rho_pfi,
    dpi = 500,
    height = 8,
    width = 14
  )
}
knitr::include_graphics(fp_p_sim_rho_pfi)
```

<img src="../figs/simulation_pfi_rho_phi.png" width="7000" />

Comparison of ZFI with varied rho and phi:

``` r
# Create/load plot
fp_p_sim_rho_zfi = "../figs/simulation_zfi_rho_phi.png"
if (!file.exists(fp_p_sim_rho_zfi)) {
  p_sim_rho_zfi <-
    agg_res_phi_rho |>
    filter(feat_imp == "ZFI", nblocks == 3) |>
    filter(t_adj <= 70) |>
    mutate(
      rhoX = factor(
        rhoX,
        levels = c("rhoX=0.2", "rhoX=0.8"),
        labels = c(expression(paste(rho[Z], "=0.2")), expression(paste(rho[Z], "=0.8")))
      ),
      rhoDelta = factor(
        rhoDelta,
        levels = c("rhoDelta=0.2", "rhoDelta=0.8"),
        labels = c(expression(paste(rho[delta], "=0.2")), expression(paste(rho[delta], "=0.8")))
      ),
      phiX = factor(
        phiX,
        levels = c("phiX=0.2", "phiX=0.8"),
        labels = c(expression(paste(phi[Z], "=0.2")), expression(paste(phi[Z], "=0.8")))
      ),
      phiDelta = factor(
        phiDelta,
        levels = c("phiDelta=0.2", "phiDelta=0.8"),
        labels = c(expression(paste(phi[delta], "=0.2")), expression(paste(phi[delta], "=0.8")))
      ),
      variable = factor(
        variable,
        levels = c("x1", "x2"),
        labels = c("Z1", "Z2")
        #labels = c(expression(Z[1]), expression(Z[2]))
      )
    ) |>
    ggplot(aes(
      x = t_adj,
      y = m_importance,
      group = interaction(variable, as.factor(simGroup)),
      color = variable
    ),
    labeller = label_parsed) +
    geom_line(alpha = 0.75) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               alpha = 0.25) +
    facet_grid(phiX + phiDelta ~ rhoX + rhoDelta, labeller = label_parsed) +
    scale_color_manual(values = cbPalette[c(7, 6)]) +
    labs(
      x = "Time",
      y = "Feature importance",
      color = "Variable",
      title = "Zeroed Feature Importances (nblocks=3)"
    ) +
    theme_bw(base_size = 18)
  ggsave(
    filename = fp_p_sim_rho_zfi,
    plot = p_sim_rho_zfi,
    dpi = 500,
    height = 8,
    width = 14
  )
}
knitr::include_graphics(fp_p_sim_rho_zfi)
```

<img src="../figs/simulation_zfi_rho_phi.png" width="7000" />
