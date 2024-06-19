# This script trains an ESN and FI to the hsw counterfactual data




library(dplyr)
library(lubridate)
library(tidyr)
library(tidync)
library(listenr)
library(stringr)
library(foreach)
library(doParallel)


# m parameter
m= 3

# number of ESN ensembles
n_ensembles = 10


# number of month blocks for FI
nblocks = 3





cl <- makeCluster(5)
registerDoParallel(cl)

fp = paste0("/gpfs/cldera/data/HSW/post/release_011423/counter_factual_latlon/")

strat_h2_file <- list.files(paste0(fp), pattern = 'T050.h2')
surf_h2_file <- list.files(paste0(fp), pattern = 'T1000.h2')

T050_raw = tidync(x = paste0(fp, strat_h2_file))
T1000_raw = tidync(x = paste0(fp, surf_h2_file))

#convert nc to tibble
T050 <-
  T050_raw %>%
  hyper_tibble() %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%


T1000 <-
  T1000_raw %>%
  hyper_tibble() %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%


hsw_all <- 
  full_join(T1000, T050,
            by = c("lon", "lat", "time")) %>%
  dplyr::select(lon, lat, time, everything()) 


# compute climatologies
hsw_clmtg_spatial_means_and_sds1 <-
  hsw_all %>%
  group_by(lon, lat) %>%
  summarise(
    T050_mean = mean(T050, na.rm = TRUE),
    T1000_mean = mean(T1000, na.rm = TRUE),
    T050_sd = sd(T050, na.rm = TRUE),
    T1000_sd = sd(T1000, na.rm = TRUE),
    .groups = "drop"
  )

# compute anomalies
hsw <-
  left_join(hsw_all, hsw_clmtg_spatial_means_and_sds1, 
            by = c("lon", "lat")) %>%
  mutate(
    T050_stndzd = (T050 - T050_mean) / (T050_sd),
    T1000_stndzd = (T1000 - T1000_mean) / (T1000_sd)
  ) %>%
  filter(time %% 10 == 0) #subsetting data for speed



hsw_train <- hsw 
hsw_train$AOD_stndzd <- rnorm(nrow(hsw_train))


# get order of lat/lon for weights
a <- hsw_train %>%
  select(time, lon, lat, T1000_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = time,
    values_from = T1000_stndzd
  ) 

lat_flat <- a$lat
weights <- cos(pi/180*lat_flat)

# training matrices
train_mat_AOD <-
  hsw_train %>%
  select(time, lon, lat, AOD_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = time,
    values_from = AOD_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_T050 <-
  hsw_train %>%
  select(time, lon, lat, T050_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = time,
    values_from = T050_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_T1000 <-
  hsw_train %>%
  select(time, lon, lat, T1000_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = time,
    values_from = T1000_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

## EOFs

n_eofs <- 5
eofs_AOD = compute_eofs(Ztrain = train_mat_AOD, n_eofs = n_eofs)
eofs_T050 = compute_eofs(Ztrain = train_mat_T050, n_eofs = n_eofs)
eofs_T1000 = compute_eofs(Ztrain = train_mat_T1000, n_eofs = n_eofs)

# Model


x_T050 = cbind(eofs_AOD$train, eofs_T050$train)
y_T050 = eofs_T050$train

t = rownames(x_T050)

# so you don't have to change for compute_fi
y_train_mat_T050 <- train_mat_T050#
phi_train_T050 <- eofs_T050$phi 

x_T1000 = cbind(eofs_AOD$train, eofs_T1000$train)
y_T1000 = eofs_T1000$train

# so you don't have to change for compute_fi 
y_train_mat_T1000 <- train_mat_T1000
phi_train_T1000 <- eofs_T1000$phi 

# train ESNs
esn_T050 <-
  fit_Eesn(
    x = x_T050,
    y = y_T050,
    t = as.character(t),
    tau = 1,
    m = m,
    tau_emb = 1,
    nh = 200,
    W_pi=0.5,
    nu=0.1,
    add_quad = TRUE,
    internal_scaling = "joint",
    seed = 20230413,
    reg_par = 50,
    n_ensembles = n_ensembles
  )

esn_T1000 <-
  fit_Eesn(
    x = x_T1000,
    y = y_T1000,
    t = as.character(t),
    tau = 1,
    m = m,
    tau_emb = 1,
    nh = 50,
    U_pi =0.5,
    W_pi=0.5,
    nu=0.1,
    add_quad = TRUE,
    internal_scaling = "joint",
    seed = 20230413,
    reg_par = 50,
    n_ensembles = n_ensembles
  )




# Compute ZFI
zfi_wo_retrain_spatial_T050 <-
  compute_fi(
    model = esn_T050,
    type="zfi",
    var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs)),
    y_spatial = y_train_mat_T050,
    phi = phi_train_T050,
    seed = 20230411,
    blockSize = nblocks,
    return_adj_preds=TRUE,
    weights=weights
  )

zfi_wo_retrain_spatial_T1000 <-
  compute_fi(
    model = esn_T1000,
    type="zfi",
    var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs)), 
    y_spatial = y_train_mat_T1000,
    phi = phi_train_T1000,
    seed = 20230411,
    blockSize = nblocks,
    return_adj_preds=TRUE,
    weights=weights
  )



# format results
zfi_res_T050 <-
zfi_wo_retrain_spatial_T050 %>%
  tibble::remove_rownames() %>%
  select(t_adj,t_forecasted,vars_adj,rmses_obs,rmses_adj,fi) %>%
  mutate(
    feat_imp = "ZFI",
    variable = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AOD","T050")
  )


zfi_res_T1000 <-
  zfi_wo_retrain_spatial_T1000 %>%
  tibble::remove_rownames() %>%
  select(t_adj,t_forecasted,vars_adj,rmses_obs,rmses_adj,fi) %>%
  mutate(
    feat_imp = "ZFI",
    variable = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AOD","T1000")
  ) 


# Feature Importance Results 

full_zfi_wo_retrain_spatial <- rbind(zfi_res_T050 %>% mutate(model="T050"),zfi_res_T1000 %>% mutate(model="T1000"))
fi_res <- full_zfi_wo_retrain_spatial 



write.csv(fi_res,file=paste0("../output/gmd_paper/107_HSW_cf.csv"),row.names=FALSE)


