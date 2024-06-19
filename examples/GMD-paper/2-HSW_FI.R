# This script trains an ESN and FI to the hsw data




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

# coutnerfactual
fp_cf = paste0("/gpfs/cldera/data/HSW/post/release_011423/counter_factual_latlon/")

strat_h2_file_cf <- list.files(paste0(fp_cf), pattern = 'T050.h2')
surf_h2_file_cf <- list.files(paste0(fp_cf), pattern = 'T1000.h2')

T050_raw_cf = tidync(x = paste0(fp_cf, strat_h2_file_cf))
T1000_raw_cf = tidync(x = paste0(fp_cf, surf_h2_file_cf))


# extract coutnerfactual nc files as tibble
T050_cf <-
  T050_raw_cf %>%
  hyper_tibble() %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon)) 


T1000_cf <-
  T1000_raw_cf %>%
  hyper_tibble() %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon)) 


hsw_all_cf <- 
  full_join(T1000_cf, T050_cf,
            by = c("lon", "lat", "time")) %>%
  dplyr::select(lon, lat, time, everything()) 

# compute climatologies on counterfactual
hsw_clmtg_spatial_means_and_sds1 <-
  hsw_all_cf %>%
  group_by(lon, lat) %>%
  summarise(
    T050_mean = mean(T050, na.rm = TRUE),
    T1000_mean = mean(T1000, na.rm = TRUE),
    T050_sd = sd(T050, na.rm = TRUE),
    T1000_sd = sd(T1000, na.rm = TRUE),
    .groups = "drop"
  )






cl <- makeCluster(5)
registerDoParallel(cl)

fp = paste0("/gpfs/cldera/data/HSW/post/release_011423/ens_members_latlon/")
sub_folder_list = list(paste0("ens0",1:5,"/"))

out <- foreach(i=1:5,.packages=c('listenr','dplyr','lubridate','tidyr','tidync','stringr')) %dopar% {
  

  sub_folder1 <- sub_folder_list[[1]][i]
  AOD_h2_file <- list.files(paste0(fp, sub_folder1), pattern = 'AOD.h2')
  strat_h2_file <- list.files(paste0(fp, sub_folder1), pattern = 'T050.h2')
  surf_h2_file <- list.files(paste0(fp, sub_folder1), pattern = 'T1000.h2')
  
  AOD_raw <- tidync(x = paste0(fp, sub_folder1, AOD_h2_file))
  T050_raw = tidync(x = paste0(fp, sub_folder1, strat_h2_file))
  T1000_raw = tidync(x = paste0(fp, sub_folder1, surf_h2_file))
  
# convert nc files to tibble
  AOD <- 
    AOD_raw %>% 
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  
  T050 <-
    T050_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  
  T1000 <-
    T1000_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  
  hsw_all <- 
    full_join(AOD, T050,
              by = c("lon", "lat", "time")) %>%
    full_join(T1000,
              by = c("lon", "lat", "time")) %>%
    dplyr::select(lon, lat, time, everything()) 
  
  hsw_clmtg_spatial_means_and_sds2 <-
    hsw_all %>%
    group_by(lon, lat) %>%
    summarise(
      AOD_mean = mean(AOD, na.rm = TRUE),
      AOD_sd = sd(AOD, na.rm = TRUE),
      .groups = "drop"
    )
  
  # compute anomalies 
  hsw <-
    left_join(hsw_all, left_join(hsw_clmtg_spatial_means_and_sds1, hsw_clmtg_spatial_means_and_sds2,by = c("lon", "lat")),
              by = c("lon", "lat")) %>%
    mutate(
      AOD_stndzd = AOD/AOD_sd,#(AOD - AOD_mean) / (AOD_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      T1000_stndzd = (T1000 - T1000_mean) / (T1000_sd)
    ) %>%
    filter(time %% 10 == 0) 
  
  
  hsw_train <- hsw 
  
  
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
  
  # create training matrices
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
  
  n_eofs <- 20
  eofs_AOD = compute_eofs(Ztrain = train_mat_AOD, n_eofs = n_eofs)
  eofs_T050 = compute_eofs(Ztrain = train_mat_T050, n_eofs = n_eofs)
  eofs_T1000 = compute_eofs(Ztrain = train_mat_T1000, n_eofs = n_eofs)
  
  # Model matrices
  
  
  x_T050 = cbind(eofs_AOD$train, eofs_T050$train)
  y_T050 = eofs_T050$train
  
  t = rownames(x_T050)
  
  # for convenience
  y_train_mat_T050 <- train_mat_T050#
  phi_train_T050 <- eofs_T050$phi 
  
  x_T1000 = cbind(eofs_AOD$train, eofs_T1000$train)
  y_T1000 = eofs_T1000$train
  
  y_train_mat_T1000 <- train_mat_T1000
  phi_train_T1000 <- eofs_T1000$phi 

  # train ESN
  #T050
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
  
  #T1000
  esn_T1000 <-
    fit_Eesn(
      x = x_T1000,
      y = y_T1000,
      t = as.character(t),
      tau = 1,
      m = m,
      tau_emb = 1,
      nh = 100,
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
  

  #format results
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
  fi_res <- full_zfi_wo_retrain_spatial %>% 
    mutate(rep=i)
  
  
  
  #--------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------
  
  #for spatial plots
  
  hsw_locations <-
    hsw_train %>%
    select(time, lon, lat, T050) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = time,
      values_from = T050
    ) %>%
    mutate(loc_id = 1:n()) %>%
    select(loc_id, lon, lat) %>%
    left_join(hsw_train,by=c('lat','lon'))%>%
    mutate(time=as.character(time))
  
  
  #Compute ESN predictions on spatial scale:
  esn_preds_spatial_T050_list <- lapply(1:length(esn_T050), function(x) predict_esn(model = esn_T050[[x]], phi = phi_train_T050)$preds_ins)
  esn_preds_spatial_T1000_list <- lapply(1:length(esn_T1000), function(x) predict_esn(model = esn_T1000[[x]], phi = phi_train_T1000)$preds_ins)
  
  esn_preds_spatial_T050 = Reduce("+", esn_preds_spatial_T050_list) / length(esn_preds_spatial_T050_list)
  esn_preds_spatial_T1000 = Reduce("+", esn_preds_spatial_T1000_list) / length(esn_preds_spatial_T1000_list)
  
  #Convert spatial predictions to long data frame:
  esn_preds_spatial_df_T050 <-
    data.frame(esn_preds_spatial_T050) %>%
    tibble::rownames_to_column(var = "time") %>%
    pivot_longer(
      cols = -time,
      names_to = "loc_id",
      values_to = "pred"
    ) %>%
    mutate(loc_id = as.integer(str_remove(loc_id, "X")),
           time=(time))
  
  esn_preds_spatial_df_T1000 <-
    data.frame(esn_preds_spatial_T1000) %>%
    tibble::rownames_to_column(var = "time") %>%
    pivot_longer(
      cols = -time,
      names_to = "loc_id",
      values_to = "pred"
    ) %>%
    mutate(loc_id = as.integer(str_remove(loc_id, "X")),
           time = (time))
  
  
  # format results
  zfi_spatial_df_T050 <-
    zfi_wo_retrain_spatial_T050 %>%
    rename(time = t_adj) %>%
    select(-c(t_forecasted,rmses_obs,rmses_adj,fi)) %>%
    pivot_longer(cols = -c(vars_adj, time), names_to = "loc_id", values_to = "pred_zeroed") %>%
    mutate(loc_id = as.integer(str_remove(loc_id, "X"))) %>%
    mutate(vars_adj = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AOD","T050"),
           rep = i
    )
  

  
  zfi_spatial_df_T1000 <-
    zfi_wo_retrain_spatial_T1000 %>%
    rename(time = t_adj) %>%
    select(-c(t_forecasted,rmses_obs,rmses_adj,fi)) %>%
    pivot_longer(cols = -c(vars_adj, time), names_to = "loc_id", values_to = "pred_zeroed") %>%
    mutate(loc_id = as.integer(str_remove(loc_id, "X"))) %>%
    mutate(vars_adj = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AOD",ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "T1000",ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"T1000","T1000"))),
           rep = i
    )
  
  
  preds_plus_T050 <-
    zfi_spatial_df_T050 %>%
    left_join(hsw_locations, by = c('time', 'loc_id'))  %>%
    left_join(esn_preds_spatial_df_T050, by = c('time', 'loc_id')) %>%
    mutate(weight = sqrt(cos(pi/180*lat))) %>%
    group_by(vars_adj, time) %>%
    mutate(total_weight = sum(weight, na.rm = T)) %>%
    ungroup() %>%
    select(vars_adj, time, loc_id, lon, lat, weight,
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
      zfi_contr_wgt = sqrt(weight * (pred_minus_obs_zeroed^2)) - sqrt(weight * (pred_minus_obs^2)),
      rep = i
    ) %>%
    group_by(lat,vars_adj,time,rep) %>%
    summarize(m_zfi = mean(zfi_contr_wgt),s_zfi=sd(zfi_contr_wgt))%>% #taking mean across ESN ensembles
    mutate(model="T050")
  
  preds_plus_T1000 <-
    zfi_spatial_df_T1000 %>%
    left_join(hsw_locations, by = c('time', 'loc_id'))  %>%
    left_join(esn_preds_spatial_df_T1000, by = c('time', 'loc_id')) %>%
    mutate(weight = sqrt(cos(pi/180*lat))) %>%
    group_by(vars_adj, time) %>%
    mutate(total_weight = sum(weight, na.rm = T)) %>%
    ungroup() %>%
    select(vars_adj, time, loc_id, lon, lat, weight,
           total_weight, T1000_stndzd, pred, pred_zeroed) %>%
    mutate(
      pred_minus_obs = pred - T1000_stndzd,
      pred_minus_obs_zeroed = pred_zeroed - T1000_stndzd
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
      zfi_contr_wgt = sqrt(weight * (pred_minus_obs_zeroed^2)) - sqrt(weight * (pred_minus_obs^2)),
      rep = i
    ) %>%
    group_by(lat,vars_adj,time,rep) %>%
    summarize(m_zfi = mean(zfi_contr_wgt),s_zfi=sd(zfi_contr_wgt))%>% #taking mean across ESN ensembles
    mutate(model="T1000")
  
  preds_plus <-  rbind(preds_plus_T050,preds_plus_T1000)
  
  
  list(fi_res=fi_res,preds_plus=preds_plus)
  #fi_res
}

stopCluster(cl)


fi <- do.call("rbind",lapply(1:length(out),function(x) out[[x]]$fi_res))
preds_plus <- do.call('rbind',lapply(1:length(out),function(x) out[[x]]$preds_plus))

write.csv(fi,file=paste0("../output/gmd_paper/106_HSW_impact.csv"),row.names=FALSE)
write.csv(preds_plus,file=paste0("../output/gmd_paper/106-hsw_fullvar_preds_plus.csv"),row.names=FALSE)


