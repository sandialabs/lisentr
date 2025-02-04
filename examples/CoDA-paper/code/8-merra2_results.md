MERRA2: FI Results
================
<br> Updated on April 17, 2024

**Overview**: This script contains MERRA-2 results for the CoDA paper.

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
library(lubridate)
library(purrr)
library(stringr)
library(tidyr)
```

Load results:

``` r
fi_eesn1 <-
  map_df(
    .x = 1:6,
    .f = function(x) {
      read.csv(paste0("../output/merra2_FI_eesn1_blk", x, ".csv"))
    }
  ) |>
  mutate(params = "Set 2", nh = 1000, reg_par = 5)

fi_eesn2 <-
  map_df(
    .x = 1:6,
    .f = function(x) {
      read.csv(paste0("../output/merra2_FI_eesn2_blk", x, ".csv"))
    }
  ) |>
  mutate(params = "Set 1", nh = 50, reg_par = 0.1)

fi_eesn3 <-
  map_df(
    .x = 1:6,
    .f = function(x) {
      read.csv(paste0("../output/merra2_FI_eesn3_blk", x, ".csv"))
    }
  ) |>
  mutate(params = "Set 3", nh = 2000, reg_par = 50)
```

Join results:

``` r
fi <-
  bind_rows(fi_eesn1, fi_eesn2, fi_eesn3) |>
  mutate(
    t_adj = as_date(t_adj),
    t_forecasted = as_date(t_forecasted),
    variable = factor(
      variable,
      levels = c("aod", "temp_strat"),
      labels = c("AOD", "Stratospheric Temperature")
    )
  ) |>
  mutate(feat_imp = paste0("st", feat_imp)) |>
  select(params, nh, reg_par, everything())
```

Specify the number of blocks to use in plots:

``` r
plot_blocks = 6
```

## Main Text Plots

Weighted feature importance (PFI and ZFI):

``` r
# Create file path
fp_p_merra_fi = "../figs/merra2_fi.png"

# Create/load plot
if (!file.exists(fp_p_merra_fi)) {
  
  fi1_blk = fi |> filter(nh == 1000, reg_par == 5, nblocks == plot_blocks)
  p_merra_fi <-
    fi1_blk |>
    ggplot() +
    geom_vline(
      xintercept = as_date("1991-06-01"),
      color = "black",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    geom_vline(
      xintercept = as_date("1982-03-29"),
      color = "gray",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    geom_line(
      mapping = aes(x = t_adj, y = fi, group = ensemble), 
      alpha = 0.5,
      color = "grey60"
    ) +
    geom_line(
      data = fi1_blk |>
        summarise(
          mean_importance = mean(fi),
          .by = c(feat_imp, nblocks, variable, t_adj)
        ),
      aes(x = t_adj, y = mean_importance),
      alpha = 1,
      linewidth = 0.5,
      color = "black"
      
    ) +
    facet_grid(feat_imp ~ variable, switch = "y") +
    labs(
      x = "Date",
      y = "Feature importance",
      color = ""
    ) +
    theme_bw(base_size = 22) +
    theme(
      strip.background = element_blank(),
      title = element_text(size = 14),
      strip.placement = "outside",
      axis.title.y = element_blank(),
      strip.text.y.left = element_text(angle = 0), 
      axis.title.x = element_text(size = 20)
    )
  
  ggsave(
    filename = fp_p_merra_fi,
    plot = p_merra_fi,
    dpi = 500,
    width = 16,
    height = 8
  )
  
}

# Load plot
knitr::include_graphics(fp_p_merra_fi)
```

<img src="../figs/merra2_fi.png" width="2560" />

## Supplemental Material Plots

Comparing number of blocks (averaged over ensembles):

``` r
# Create file path
fp_p_merra_fi_blocks_ens_ave = "../figs/merra2_fi_blocks_ens_aves.png"

# Create/load plot
if (!file.exists(fp_p_merra_fi_blocks_ens_ave)) {
  p_merra_fi_blocks_ens_ave <-
    fi |>
    filter(nh == 1000, reg_par == 5) |>
    summarise(
      mean_importance = mean(fi),
      .by = c(feat_imp, nblocks, variable, t_adj)
    ) |>
    ggplot(aes(
      x = t_adj,
      y = mean_importance,
      group = nblocks,
      color = as.factor(nblocks)
    )) +
    geom_vline(
      xintercept = as_date("1991-06-01"),
      color = "black",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_vline(
      xintercept = as_date("1982-03-29"),
      color = "gray",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_line(alpha = 0.85) +
    facet_grid(feat_imp ~ variable, switch = "y") +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1, "black")) +
    labs(
      x = "Date",
      y = "Feature importance",
      color = "Block size"
    ) +
    theme_bw(base_size = 18) +
    theme(
      strip.background = element_blank(),
      title = element_text(size = 14),
      strip.placement = "outside",
      axis.title.y = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(
      nrow = 1,
      override.aes = list(alpha = 1, linewidth = 1)
    ))
  ggsave(
    filename = fp_p_merra_fi_blocks_ens_ave,
    plot = p_merra_fi_blocks_ens_ave,
    dpi = 500,
    width = 14,
    height = 8
  )
}

# Load plot
knitr::include_graphics(fp_p_merra_fi_blocks_ens_ave)
```

<img src="../figs/merra2_fi_blocks_ens_aves.png" width="2240" />

Comparing number of blocks:

``` r
# Create file path
fp_p_merra_fi_blocks = "../figs/merra2_fi_blocks_ens.png"

# Create/load plot
if (!file.exists(fp_p_merra_fi_blocks)) {
  
  p_merra_fi_blocks <-
    fi |>
    filter(nh == 1000, reg_par == 5) |>
    ggplot(aes(
      x = t_adj,
      y = fi,
      group = ensemble,
      color = as.factor(nblocks)
    )) +
    geom_vline(
      xintercept = as_date("1991-06-01"),
      color = "black",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_vline(
      xintercept = as_date("1982-03-29"),
      color = "gray",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_line(alpha = 0.15) +
    facet_grid(nblocks ~ paste0(variable, ": ", feat_imp), switch = "y") +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1, "black")) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 0.5))) +
    labs(x = "Date",
         y = "Feature importance",
         color = "Block size") +
    theme_bw(base_size = 24) +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      title = element_text(size = 14),
      strip.placement = "outside",
      legend.title = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_blank(),
      strip.text.y.left = element_blank()
    ) +
    guides(color = guide_legend(
      nrow = 1,
      override.aes = list(alpha = 1, linewidth = 1)
    ))
  
  ggsave(
    filename = fp_p_merra_fi_blocks,
    plot = p_merra_fi_blocks,
    dpi = 500,
    width = 20,
    height = 12
  )
  
}

# Load plot
knitr::include_graphics(fp_p_merra_fi_blocks)
```

<img src="../figs/merra2_fi_blocks_ens.png" width="3200" />

Comparing hyperparameter values (averaged over ensembles):

``` r
# Create file path
fp_p_fi_params = "../figs/merra2_fi_param_comparison.png"

# Create/load plot
if (!file.exists(fp_p_fi_params)) {
  p_fi_params <-
    fi |>
    filter(nblocks == plot_blocks) |>
    summarise(
      mean_importance = mean(fi),
      .by = c(params, feat_imp, nblocks, variable, t_adj)
    ) |>
    mutate(params = factor(params)) |>
    ggplot() +
    geom_vline(
      xintercept = as_date("1991-06-01"),
      color = "black",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    geom_vline(
      xintercept = as_date("1982-03-29"),
      color = "gray",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    geom_line(
      aes(
        x = t_adj,
        y = mean_importance,
        color = params,
        group = params
      ),
      alpha = 0.9,
      linewidth = 0.75
    ) +
    facet_grid(feat_imp ~ variable, switch = "y") +
    labs(x = "Date",
         y = "Feature importance",
         color = "Parameters") +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1)) +
    theme_bw(base_size = 18) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.title.y = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      legend.title = element_text(size = 16)
    ) +
    guides(color =
             guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
  
  ggsave(
    filename = fp_p_fi_params,
    plot = p_fi_params,
    dpi = 500,
    width = 14,
    height = 8
  )
  
}

# Load plot
knitr::include_graphics(fp_p_fi_params)
```

<img src="../figs/merra2_fi_param_comparison.png" width="2240" />

Comparing hyperparameter values:

``` r
# Create file path
fp_p_fi_params_var = "../figs/merra2_fi_param_var_comparison.png"

# Create/load plot
if (!file.exists(fp_p_fi_params_var)) {
  p_fi_params_var <-
    fi |>
    mutate(params = factor(params)) |>
    filter(nblocks == plot_blocks) |>
    ggplot() +
    geom_vline(
      xintercept = as_date("1991-06-01"),
      color = "black",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    geom_vline(
      xintercept = as_date("1982-03-29"),
      color = "gray",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    geom_line(
      mapping = aes(
        x = t_adj,
        y = fi,
        group = ensemble,
        color = params
      ),
      alpha = 0.5,
    ) +
    facet_grid(paste0(feat_imp, "\n ", params) ~ variable, switch = "y") +
    labs(x = "Date",
         y = "Feature importance",
         color = "Parameters") +
    scale_color_manual(values = c(wesanderson::wes_palettes$Zissou1)) +
    theme_bw(base_size = 18) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.title.y = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      legend.title = element_text(size = 16)
    ) + 
    guides(color =
             guide_legend(override.aes = list(alpha = 1, linewidth = 0.75)))
  ggsave(
    filename = fp_p_fi_params_var,
    plot = p_fi_params_var,
    dpi = 500,
    width = 14,
    height = 8
  )
  
}

# Load plot
knitr::include_graphics(fp_p_fi_params_var)
```

<img src="../figs/merra2_fi_param_var_comparison.png" width="2240" />

## Plot Not Included in Manuscript

Comparison of RMSEs before and after permuted/zeroed:

``` r
# Create file path
fp_p_merra_fi_rmses = "../figs/merra2_fi_rmses.png"

# Create/load plot
if (!file.exists(fp_p_merra_fi_rmses)) {

  fi_rmses <- 
    fi |>
    filter(nh == 1000, reg_par == 5) |>
    select(ensemble, feat_imp, nblocks, variable, t_adj,
           rmses_obs, rmses_adj) |>
    pivot_longer(cols = rmses_obs:rmses_adj,
                 names_to = "rmse_type",
                 values_to = "rmse") |>
    mutate(rmse_type = str_remove(rmse_type, "rmses_")) |>
    filter(nblocks == plot_blocks)
    
  p_merra_fi_rmses <-
    fi_rmses |>
    ggplot(aes(x = t_adj, color = rmse_type)) +
    geom_vline(
      xintercept = as_date("1991-06-01"),
      color = "black",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    geom_vline(
      xintercept = as_date("1982-03-29"),
      color = "gray",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    geom_line(aes(y = rmse, group = factor(ensemble):factor(rmse_type)), 
              alpha = 0.5, linewidth = 0.1) +
    geom_line(
      data = fi_rmses  |>
      summarise(
        mean_rmse = mean(rmse),
        .by = c(feat_imp, variable, t_adj, rmse_type)
      ),
      aes(y = mean_rmse),
      alpha = 1,
      linewidth = 0.75
    ) +
    facet_grid(feat_imp ~ variable) +
    scale_color_manual(
      values = c(wesanderson::wes_palettes$Zissou1[5], "black")
    ) +
    labs(x = "Date",
         y = "RMSE (weighted)",
         color = "") +
    theme_bw(base_size = 18) +
    theme(
      strip.background = element_blank(),
      title = element_text(size = 14),
      strip.text.y.right = element_text(angle = 0),
      legend.position = "bottom"
    )
  
  ggsave(
    filename = fp_p_merra_fi_rmses,
    plot = p_merra_fi_rmses,
    dpi = 500,
    width = 14,
    height = 8
  )
  
}

# Load plot
knitr::include_graphics(fp_p_merra_fi_rmses)
```

<img src="../figs/merra2_fi_rmses.png" width="2240" />
