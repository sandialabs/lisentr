# this script trains a ESN and FI on MERRA2

library(dplyr)
library(doParallel)
library(foreach)
library(ggplot2)
library(listenr)
library(lubridate)
library(purrr)
library(tidyr)
library(stringr)


  
fp = "/projects/cldera/obs_thrust/MERRA-2/"


## Raw data

TOTEXTTAU_raw = read.csv(paste0(fp, "tavgM_2d_aer_Nx/TOTEXTTAUall_360x180.csv"))
T050_raw = read.csv(paste0(fp, "instM_3d_asm_Np/MERRA2_3dasm_temperature_50mb_360x180.csv"))
TREFHT_raw = read.csv(paste0(fp, "tavg1_2d_slv_Nx_daily/MERRA2_tavg1_2d_slv_monthly_Nx_360x180.csv"))
rad_raw <- read.csv(paste0(fp,"tavgM_2d_rad_Nx/Rad_360x180.csv"))

# number of ensembles 
n_ensembles = 10
m=3 

## Cleaning

TOTEXTTAU <- 
  TOTEXTTAU_raw %>% 
  mutate(date = as_date(date))

rad <- 
  rad_raw %>% 
  mutate(date = as_date(date))

T050 <-
  T050_raw %>% 
  rename("T050" = `T`) %>%
  mutate(date = as_date(date)) %>%
  mutate(day = day(date))

TREFHT <-
  TREFHT_raw %>% 
  rename("TREFHT" = `T`) %>%
  mutate(date = as_date(date)) %>%
  mutate(day = day(date),
         month = month(date),
         year = year(date))



# combine data
merra2 <- 
  full_join(TOTEXTTAU, T050 %>% select(-c(month,year,day)), by = c("lon", "lat", "date")) %>%
  full_join(rad %>% filter(date <= as_date("1998-12-31")),by= c("lon", "lat", "date")) %>%
  full_join(TREFHT %>% select(-c(month,year,day)),by= c("lon", "lat", "date")) %>%
  mutate(year=year(date),month=month(date),data='train') %>%
  filter(year >= 1991 )




# compute climatologies
merra2_clmtg_spatial_means_and_sds <-
  merra2 %>%
  group_by(month, lon, lat) %>%
  summarise(
    TOTEXTTAU_mean = mean(TOTEXTTAU, na.rm = TRUE),
    T050_mean = mean(T050, na.rm = TRUE),
    TREFHT_mean = mean(TREFHT, na.rm = TRUE),
    SWGDNCLR_mean = mean(SWGDNCLR, na.rm = TRUE),
    LWTUP_mean = mean(LWTUP, na.rm = TRUE),
    SWGDNCLR_mean = mean(SWGDNCLR,na.rm=TRUE),
    TOTEXTTAU_sd = sd(TOTEXTTAU, na.rm = TRUE),
    T050_sd = sd(T050, na.rm = TRUE),
    TREFHT_sd = sd(TREFHT, na.rm = TRUE),
    SWGDNCLR_sd = sd(SWGDNCLR, na.rm = TRUE),
    LWTUP_sd = sd(LWTUP, na.rm = TRUE),
    SWGDNCLR_sd = sd(SWGDNCLR,na.rm=TRUE),
    .groups = "drop"
  )




  
# compute anomalies
merra2_esn <-
  left_join(merra2, merra2_clmtg_spatial_means_and_sds, by = c("month", "lon", "lat")) %>%
  filter(abs(lat) <= 67.5) %>% # remove AEROD_v NAs
  mutate(
    TOTEXTTAU_stndzd = (TOTEXTTAU - TOTEXTTAU_mean) / (TOTEXTTAU_sd),
    LWTUP_stndzd = (LWTUP - LWTUP_mean) / (LWTUP_sd),
    SWGDNCLR_stndzd = (SWGDNCLR - SWGDNCLR_mean),# / (SWGDNCLR_sd),
    SWGDNCLR_stndzd = (SWGDNCLR - SWGDNCLR_mean) / (SWGDNCLR_sd),
    T050_stndzd = (T050 - T050_mean) / (T050_sd),
    TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd)
  )


## ESN Matrices

# get order of lat/lon for weights
a <- merra2_esn %>%
  select(date, lon, lat, TREFHT_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = TREFHT_stndzd
  ) 

lat_flat <- a$lat
weights <- cos(pi/180*lat_flat)

train_mat_TOTEXTTAU <-
  merra2_esn %>%
  select(date, lon, lat, TOTEXTTAU_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = TOTEXTTAU_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()



train_mat_LWTUP <-
  merra2_esn %>%
  select(date, lon, lat, LWTUP_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = LWTUP_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_SWGDNCLR <-
  merra2_esn %>%
  select(date, lon, lat, SWGDNCLR_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = SWGDNCLR_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_SWGDNCLR <-
  merra2_esn %>%
  select(date, lon, lat, SWGDNCLR_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = SWGDNCLR_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_T050 <-
  merra2_esn %>%
  select(date, lon, lat, T050_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = T050_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_TREFHT <-
  merra2_esn %>%
  select(date, lon, lat, TREFHT_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = TREFHT_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

## EOFs

n_eofs <- 20
eofs_TOTEXTTAU = compute_eofs(Ztrain = train_mat_TOTEXTTAU, n_eofs = n_eofs)
eofs_SWGDNCLR = compute_eofs(Ztrain = train_mat_SWGDNCLR, n_eofs = n_eofs)
eofs_LWTUP = compute_eofs(Ztrain = train_mat_LWTUP, n_eofs = n_eofs)
eofs_T050 = compute_eofs(Ztrain = train_mat_T050, n_eofs = n_eofs)
eofs_TREFHT = compute_eofs(Ztrain = train_mat_TREFHT, n_eofs = n_eofs)



# Model matrices
  
x_T050 = cbind(eofs_TOTEXTTAU$train, eofs_LWTUP$train, eofs_T050$train)
y_T050 = eofs_T050$train

x_TREFHT = cbind(eofs_TOTEXTTAU$train, eofs_SWGDNCLR$train, eofs_TREFHT$train)
y_TREFHT = eofs_TREFHT$train

# for convenience
y_train_mat_T050 <- train_mat_T050 
phi_train_T050 <- eofs_T050$phi 

y_train_mat_TREFHT <- train_mat_TREFHT
phi_train_TREFHT <- eofs_TREFHT$phi 

t = rownames(x_T050)


# train ESNs
# T050
esn_T050 <-
  fit_Eesn(
    x = x_T050,
    y = y_T050,
    t = as.character(t),
    tau = 1,
    m = 3,
    tau_emb = 1,
    U_pi=0.5,
    W_pi=0.5,
    nh = 200,
    add_quad = TRUE,
    internal_scaling = "joint",
    seed = 20230413,
    reg_par = 50,
    nu=.1,
    n_ensembles=n_ensembles
  )

#TREFHT
esn_TREFHT <-
  fit_Eesn(
    x = x_TREFHT,
    y = y_TREFHT,
    t = as.character(t),
    tau = 1,
    m = 3,
    tau_emb = 1,
    nh = 200,
    add_quad = TRUE,
    internal_scaling = "joint",
    seed = 20230413,
    reg_par = 5,
    nu=.1,
    n_ensembles=n_ensembles
  )



# Compute ZFI
nblocks = 4


#T050
zfi_T050 <-
  compute_fi(
    model = esn_T050,
    type="zfi",
    var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs),(2*n_eofs + 1):(3 * n_eofs)), 
    y_spatial = y_train_mat_T050,
    phi = phi_train_T050,
    seed = 20230411,
    blockSize = nblocks,
    return_adj_preds=TRUE,
    weights=weights
  )

#TREFHT
zfi_TREFHT <-
  compute_fi(
    model = esn_TREFHT,
    type="zfi",
    var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs),(2*n_eofs + 1):(3 * n_eofs)),
    y_spatial = y_train_mat_TREFHT,
    phi = phi_train_TREFHT,
    seed = 20230411,
    blockSize = nblocks,
    return_adj_preds=TRUE,
    weights=weights
  )



# format results
zfi_res_T050 <-
  zfi_T050 %>%
  tibble::remove_rownames() %>%
  select(t_adj,t_forecasted,vars_adj,rmses_obs,rmses_adj,fi) %>%
  mutate(
    variable = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "TOTEXTTAU",ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "LWTUP",ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"T050","TREFHT"))), 
    t_forecasted = as_date(t_forecasted),
    t_adj = as_date(t_adj)
  ) %>%
  mutate(month_diff = interval(t_adj, t_forecasted) %/% months(1)) %>%
  filter(month_diff == 1) 


zfi_res_TREFHT <-
  zfi_TREFHT %>%  
  tibble::remove_rownames()  %>%
  select(t_adj,t_forecasted,vars_adj,rmses_obs,rmses_adj,fi) %>%
  mutate(
    variable = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "TOTEXTTAU",ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "SWGDNCLR",ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"TREFHT","TREFHT"))), 
    t_forecasted = as_date(t_forecasted),
    t_adj = as_date(t_adj)
  ) %>%
  mutate(month_diff = interval(t_adj, t_forecasted) %/% months(1)) %>%
  filter(month_diff == 1) 

full_zfi <- rbind(zfi_res_T050 %>% mutate(model="T050"),zfi_res_TREFHT %>% mutate(model="TREFHT"))

write.csv(full_zfi,"../output/gmd_paper/88-MERRA2_360x180_20pcs_swgndclr.csv",row.names=FALSE)



# Spatial FI plots


#for spatial plots

locations <- 
  merra2_esn %>%
  select(date, lon, lat, T050) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = T050
  ) %>%
  mutate(loc_id = 1:n()) %>%
  select(loc_id, lon, lat)

merra2_locations <- merra2_esn %>%
  left_join(locations,by=c('lat','lon'))


#
#Compute ESN predictions on spatial scale:
esn_preds_spatial_T050_list <- lapply(1:length(esn_T050), function(x) predict_esn(model = esn_T050[[x]], phi = phi_train_T050)$preds_ins)
esn_preds_spatial_TREFHT_list <- lapply(1:length(esn_TREFHT), function(x) predict_esn(model = esn_TREFHT[[x]], phi = phi_train_TREFHT)$preds_ins)

esn_preds_spatial_T050 = Reduce("+", esn_preds_spatial_T050_list) / length(esn_preds_spatial_T050_list)
esn_preds_spatial_TREFHT = Reduce("+", esn_preds_spatial_TREFHT_list) / length(esn_preds_spatial_TREFHT_list)

#Convert spatial predictions to long data frame:
esn_preds_spatial_df_T050 <-
  data.frame(esn_preds_spatial_T050) %>%
  tibble::rownames_to_column(var = "date") %>%
  pivot_longer(
    cols = -date,
    names_to = "loc_id",
    values_to = "pred"
  ) %>%
  mutate(loc_id = as.integer(str_remove(loc_id, "X")),
         date=as_date(date)) 

esn_preds_spatial_df_TREFHT <-
  data.frame(esn_preds_spatial_TREFHT) %>%
  tibble::rownames_to_column(var = "date") %>%
  pivot_longer(
    cols = -date,
    names_to = "loc_id",
    values_to = "pred"
  ) %>%
  mutate(loc_id = as.integer(str_remove(loc_id, "X")),
         date = as_date(date)) 

# format results

zfi_spatial_df_T050 <-
  zfi_T050 %>%
  rename(date = t_adj) %>%
  select(-c(t_forecasted,rmses_obs,rmses_adj,fi)) %>%
  pivot_longer(cols = -c(vars_adj, date), names_to = "loc_id", values_to = "pred_zeroed") %>%
  mutate(loc_id = as.integer(str_remove(loc_id, "X"))) %>%
  mutate(vars_adj = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "TOTEXTTAU",ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "LWTUP",ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"T050","TREFHT"))))

zfi_spatial_df_T050$date <- as_date(zfi_spatial_df_T050$date)


zfi_spatial_df_TREFHT <-
  zfi_TREFHT %>%
  rename(date = t_adj) %>%
  select(-c(t_forecasted,rmses_obs,rmses_adj,fi)) %>%
  pivot_longer(cols = -c(vars_adj, date), names_to = "loc_id", values_to = "pred_zeroed") %>%
  mutate(loc_id = as.integer(str_remove(loc_id, "X"))) %>%
  mutate(vars_adj = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "TOTEXTTAU",ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "SWGDNCLR",ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"TREFHT","TREFHT"))))

zfi_spatial_df_TREFHT$date <- as_date(zfi_spatial_df_TREFHT$date)


preds_plus_T050 <-
  zfi_spatial_df_T050 %>%
  left_join(merra2_locations, by = c('date', 'loc_id'))  %>%
  left_join(esn_preds_spatial_df_T050, by = c('date', 'loc_id')) %>%
  mutate(weight = sqrt(cos(pi/180*lat))) %>%
  group_by(vars_adj, date) %>%
  mutate(total_weight = sum(weight, na.rm = T)) %>%
  ungroup() %>%
  select(vars_adj, date, loc_id, lon, lat, weight,
         total_weight, T050_stndzd, pred, pred_zeroed) %>%
  mutate(
    pred_minus_obs = pred - T050_stndzd,
    pred_minus_obs_zeroed = pred_zeroed - T050_stndzd
  ) %>%
  mutate(
    pred_minus_obs_wgt = weight * pred_minus_obs,
    pred_minus_obs_zeroed_wgt = weight * pred_minus_obs_zeroed,
    pred_minus_obs_abs = abs(pred_minus_obs),
    pred_minus_obs_zeroed_abs = abs(pred_minus_obs_zeroed),
    pred_minus_obs_abs_wgt = weight * abs(pred_minus_obs),
    pred_minus_obs_zeroed_abs_wgt = weight * abs(pred_minus_obs_zeroed)
  ) %>%
  mutate(
    zfi_contr = pred_minus_obs_zeroed_abs - pred_minus_obs_abs,
    zfi_contr_wgt = sqrt(weight * (pred_minus_obs_zeroed^2)) - sqrt(weight * (pred_minus_obs^2))
  ) %>%
  summarize(m_zfi_contr_wgt = mean(zfi_contr_wgt),s_zfi_contr_wgt=sd(zfi_contr_wgt),.by=c(lat,vars_adj,date))%>% #taking mean across ESN ensembles
  mutate(model="T050")

preds_plus_TREFHT <-
  zfi_spatial_df_TREFHT %>%
  left_join(merra2_locations, by = c('date', 'loc_id'))  %>%
  left_join(esn_preds_spatial_df_TREFHT, by = c('date', 'loc_id')) %>%
  mutate(weight = sqrt(cos(pi/180*lat))) %>%
  group_by(vars_adj, date) %>%
  mutate(total_weight = sum(weight, na.rm = T)) %>%
  ungroup() %>%
  select(vars_adj, date, loc_id, lon, lat, weight,
         total_weight, TREFHT_stndzd, pred, pred_zeroed) %>%
  mutate(
    pred_minus_obs = pred - TREFHT_stndzd,
    pred_minus_obs_zeroed = pred_zeroed - TREFHT_stndzd
  ) %>%
  mutate(
    pred_minus_obs_wgt = weight * pred_minus_obs,
    pred_minus_obs_zeroed_wgt = weight * pred_minus_obs_zeroed,
    pred_minus_obs_abs = abs(pred_minus_obs),
    pred_minus_obs_zeroed_abs = abs(pred_minus_obs_zeroed),
    pred_minus_obs_abs_wgt = weight * abs(pred_minus_obs),
    pred_minus_obs_zeroed_abs_wgt = weight * abs(pred_minus_obs_zeroed)
  ) %>%
  mutate(
    zfi_contr = pred_minus_obs_zeroed_abs - pred_minus_obs_abs,
    zfi_contr_wgt = sqrt(weight * (pred_minus_obs_zeroed^2)) - sqrt(weight * (pred_minus_obs^2))
  ) %>%
  summarize(m_zfi_contr_wgt = mean(zfi_contr_wgt),s_zfi_contr_wgt=sd(zfi_contr_wgt),.by=c(lat,vars_adj,date))%>% #taking mean across ESN ensembles
  mutate(model="TREFHT")


preds_plus <-  rbind(preds_plus_T050 ,preds_plus_TREFHT )

write.csv(preds_plus,"../output/gmd_paper/88-MERRA2_preds_plus_20pcs_swgndclr.csv",row.names=FALSE)












