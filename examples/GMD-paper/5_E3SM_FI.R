# This script trains an ESN and FI to the E3SM fullvar 8 year run data




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
nblocks = 4

cl <- makeCluster(5)
registerDoParallel(cl)

out <- foreach(i=1:5,.packages=c('listenr','dplyr','lubridate','tidyr','tidync','stringr')) %dopar% {
  
  fp = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,"/post/atm/180x360_aave/ts/monthly/8yr/")
  fp_cf = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,".cf/post/atm/180x360_aave/ts/monthly/8yr/")
  
  AEROD_v_h2_file <- paste0(fp, "AEROD_v_199101_199812.nc")
  T050_h2_file <- paste0(fp, "T050_199101_199812.nc")
  TREFHT_h2_file <- paste0(fp, "TREFHT_199101_199812.nc")
  FLNT_h2_file <- paste0(fp, "FLNT_199101_199812.nc")
  FSDSC_h2_file <- paste0(fp, "FSDSC_199101_199812.nc")
  
  #load nc files
  AEROD_v_raw <- tidync(AEROD_v_h2_file)
  T050_raw = tidync(T050_h2_file)
  TREFHT_raw = tidync(TREFHT_h2_file)
  FLNT_raw = tidync(FLNT_h2_file)
  FSDSC_raw = tidync(FSDSC_h2_file)
  
  #coutnerfactual for baseline
  AEROD_v_h2_file_cf <- paste0(fp_cf, "AEROD_v_199101_199812.nc")
  T050_h2_file_cf <- paste0(fp_cf, "T050_199101_199812.nc")
  TREFHT_h2_file_cf <- paste0(fp_cf, "TREFHT_199101_199812.nc")
  FLNT_h2_file_cf <- paste0(fp_cf, "FLNT_199101_199812.nc")
  FSDSC_h2_file_cf <- paste0(fp_cf, "FSDSC_199101_199812.nc")
  
  # load nc files
  AEROD_v_raw_cf <- tidync(AEROD_v_h2_file_cf)
  T050_raw_cf = tidync(T050_h2_file_cf)
  TREFHT_raw_cf = tidync(TREFHT_h2_file_cf)
  FLNT_raw_cf = tidync(FLNT_h2_file_cf)
  FSDSC_raw_cf = tidync(FSDSC_h2_file_cf)
  
  
  #convert nc files to tibble
  
  AEROD_v <- 
    AEROD_v_raw %>% 
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  T050 <-
    T050_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  
  TREFHT <-
    TREFHT_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  
  FLNT <-
    FLNT_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  
  FSDSC <-
    FSDSC_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  #combine data
  e3sm_all <- 
    full_join(AEROD_v, T050,
              by = c("lon", "lat", "time")) %>%
    full_join(TREFHT,
              by = c("lon", "lat", "time")) %>%
    full_join(FLNT,
              by = c("lon", "lat", "time")) %>%
    full_join(FSDSC,
              by = c("lon", "lat", "time")) %>%
    dplyr::select(lon, lat, time, everything()) %>%
    mutate(
      date = as_date("1991-01-31") + days(time-min(time)),
      month=month(date)
    )
  
  
  
  #counterfactuals
  AEROD_v_cf <- 
    AEROD_v_raw_cf %>% 
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) # %>%
  
  
  T050_cf <-
    T050_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  
  TREFHT_cf <-
    TREFHT_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  FLNT_cf <-
    FLNT_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  
  FSDSC_cf <-
    FSDSC_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  # combine
  e3sm_all_cf <- 
    full_join(AEROD_v_cf, T050_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(TREFHT_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(FLNT_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(FSDSC_cf,
              by = c("lon", "lat", "time")) %>%
    dplyr::select(lon, lat, time, everything()) %>%
    mutate(
      date = as_date("1991-01-31") + days(time-min(time)),
      month=month(date)
    )
  
  
  # compute climatologies
  e3sm_clmtg_spatial_means_and_sds1 <-
    e3sm_all_cf %>%
    group_by(lon, lat,month) %>%
    summarise(
      AEROD_v_mean = mean(AEROD_v, na.rm = TRUE),
      T050_mean = mean(T050, na.rm = TRUE),
      TREFHT_mean = mean(TREFHT, na.rm = TRUE),
      FLNT_mean = mean(FLNT, na.rm = TRUE),
      FSDSC_mean = mean(FSDSC, na.rm = TRUE),
      AEROD_v_sd = sd(AEROD_v, na.rm = TRUE),
      T050_sd = sd(T050, na.rm = TRUE),
      TREFHT_sd = sd(TREFHT, na.rm = TRUE),
      FLNT_sd = sd(FLNT, na.rm = TRUE),
      FSDSC_sd = sd(FSDSC, na.rm = TRUE),
      .groups = "drop"
    )
  
  # compute anomalies
  e3sm <-
    left_join(e3sm_all, e3sm_clmtg_spatial_means_and_sds1, 
              by = c("lon", "lat","month")) %>%
    mutate(
      AEROD_v_stndzd = (AEROD_v - AEROD_v_mean) / (AEROD_v_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd),
      FLNT_stndzd = (FLNT - FLNT_mean) / (FLNT_sd),
      FSDSC_stndzd = (FSDSC - FSDSC_mean) / (FSDSC_sd)
      
    )
  
  
  #-------------------
  #remove AOD NAs
  
  b1=e3sm %>% group_by(lon,lat) %>% summarize(n=sum(is.na(AEROD_v)),n2=sum(is.na(FSDSC)),n3=n+n2)
  
  complete_latlon <- b1[b1$n3==0,1:2]
  
  #---------------
  
  
  e3sm_train <- e3sm %>%
    inner_join(complete_latlon)
  
  
  # get order of lat/lon for weights
  a <- e3sm_train %>%
    select(date, lon, lat, TREFHT_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = TREFHT_stndzd
    ) 
  
  lat_flat <- a$lat
  weights <- cos(pi/180*lat_flat)
  
  # train matrices
  train_mat_AEROD_v <-
    e3sm_train %>%
    select(date, lon, lat, AEROD_v_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = AEROD_v_stndzd
    ) %>%
    select(-lon,-lat) %>%
    t()
  
  
  
  train_mat_FLNT <-
    e3sm_train %>%
    select(date, lon, lat, FLNT_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = FLNT_stndzd
    ) %>%
    select(-lon,-lat) %>%
    t()
  
  train_mat_FSDSC <-
    e3sm_train %>%
    select(date, lon, lat, FSDSC_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = FSDSC_stndzd
    ) %>%
    select(-lon,-lat) %>%
    t()
  
  train_mat_T050 <-
    e3sm_train %>%
    select(date, lon, lat, T050_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = T050_stndzd
    ) %>%
    select(-lon,-lat) %>%
    t()
  
  train_mat_TREFHT <-
    e3sm_train %>%
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
  eofs_AEROD_v = compute_eofs(Ztrain = train_mat_AEROD_v, n_eofs = n_eofs)
  eofs_FSDSC = compute_eofs(Ztrain = train_mat_FSDSC, n_eofs = n_eofs)
  eofs_FLNT = compute_eofs(Ztrain = train_mat_FLNT, n_eofs = n_eofs)
  eofs_T050 = compute_eofs(Ztrain = train_mat_T050, n_eofs = n_eofs)
  eofs_TREFHT = compute_eofs(Ztrain = train_mat_TREFHT, n_eofs = n_eofs)
  
  
  
  # Model matrices
  x_T050 = cbind(eofs_AEROD_v$train, eofs_FLNT$train, eofs_T050$train)
  y_T050 = eofs_T050$train
  
  # for convenience 
  y_train_mat_T050 <- train_mat_T050#
  phi_train_T050 <- eofs_T050$phi 
  
  x_TREFHT = cbind(eofs_AEROD_v$train, eofs_FSDSC$train, eofs_TREFHT$train)
  y_TREFHT = eofs_TREFHT$train
  
  y_train_mat_TREFHT <- train_mat_TREFHT#TREFHT #train_mat_TREFHT 
  phi_train_TREFHT <- eofs_TREFHT$phi #eofs_T050$phi
  
  
  t = rownames(x_T050)
  
  
  #train ESN
  #T050
  esn_T050 <-
    fit_Eesn(
      x = x_T050,
      y = y_T050,
      t = as.character(t),
      tau = 1,
      m = m,
      tau_emb = 1,
      nh = 50,
      U_pi=0.5,
      W_pi=0.5,
      add_quad = FALSE,
      internal_scaling = "joint",
      seed = 20230413,
      reg_par = 5,
      nu=.1,
      n_ensembles = n_ensembles
    )
  
  #TREFHT
  esn_TREFHT <-
    fit_Eesn(
      x = x_TREFHT,
      y = y_TREFHT,
      t = as.character(t),
      tau = 1,
      m = m,
      tau_emb = 1,
      nh = 200,
      add_quad = FALSE,
      internal_scaling = "joint",
      seed = 20230413,
      reg_par = 5,
      nu=.1,
      n_ensembles = n_ensembles
    )
  
  
  
  
  
  # Compute ZFI
# T050
  zfi_wo_retrain_spatial_T050 <-
    compute_fi(
      model = esn_T050,
      type="zfi",
      var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs),(2*n_eofs + 1):(3 * n_eofs)),#,(3*n_eofs + 1):(4 * n_eofs)), 
      y_spatial = y_train_mat_T050,
      phi = phi_train_T050,
      seed = 20230411,
      blockSize = nblocks,
      return_adj_preds=TRUE,
      weights=weights
    )
  
  #TREFHT
  zfi_wo_retrain_spatial_TREFHT <-
    compute_fi(
      model = esn_TREFHT,
      type="zfi",
      var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs),(2*n_eofs + 1):(3 * n_eofs)),#,(3*n_eofs + 1):(4 * n_eofs)), 
      y_spatial = y_train_mat_TREFHT,
      phi = phi_train_TREFHT,
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
      variable = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AEROD_v",ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "FLNT",ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"T050","TREFHT"))), 
      t_forecasted = as_date(t_forecasted),
      t_adj = as_date(t_adj)
    ) %>%
    mutate(month_diff = interval(floor_date(as_date(t_adj),'month'),floor_date(as_date(t_forecasted),'month')) %/% months(1)) %>%
    filter(month_diff == 1) 
  
  
  zfi_res_TREFHT <-
    zfi_wo_retrain_spatial_TREFHT %>%
    tibble::remove_rownames() %>%
    select(t_adj,t_forecasted,vars_adj,rmses_obs,rmses_adj,fi) %>%
    mutate(
      feat_imp = "ZFI",
      variable = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AEROD_v",ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "FSDSC",ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"TREFHT","TREFHT"))), 
      t_forecasted = as_date(t_forecasted),
      t_adj = as_date(t_adj)
    ) %>%
    mutate(month_diff = interval(floor_date(as_date(t_adj),'month'),floor_date(as_date(t_forecasted),'month')) %/% months(1)) %>%
    filter(month_diff == 1)   
  
  
  # Feature Importance Results 
  
  full_zfi_wo_retrain_spatial <- rbind(zfi_res_T050 %>% mutate(model="T050"),zfi_res_TREFHT %>% mutate(model="TREFHT"))
  fi_res <- full_zfi_wo_retrain_spatial %>% 
    mutate(rep=i)
  
  
  
  #--------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------
  
  #for spatial plots
  
  e3sm_locations <- 
    e3sm_train %>%
    select(date, lon, lat, T050) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = T050
    ) %>%
    mutate(loc_id = 1:n()) %>%
    select(loc_id, lon, lat) %>%
    left_join(e3sm_train,by=c('lat','lon'))
  

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
    zfi_wo_retrain_spatial_T050 %>%
    rename(date = t_adj) %>%
    select(-c(t_forecasted,rmses_obs,rmses_adj,fi)) %>%
    pivot_longer(cols = -c(vars_adj, date), names_to = "loc_id", values_to = "pred_zeroed") %>%
    mutate(loc_id = as.integer(str_remove(loc_id, "X"))) %>%
    mutate(vars_adj = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AEROD_v",ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "FLNT",ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"T050","TREFHT"))),
           rep = i
    )
  
  zfi_spatial_df_T050$date <- as_date(zfi_spatial_df_T050$date)
  
  
  zfi_spatial_df_TREFHT <-
    zfi_wo_retrain_spatial_TREFHT %>%
    rename(date = t_adj) %>%
    select(-c(t_forecasted,rmses_obs,rmses_adj,fi)) %>%
    pivot_longer(cols = -c(vars_adj, date), names_to = "loc_id", values_to = "pred_zeroed") %>%
    mutate(loc_id = as.integer(str_remove(loc_id, "X"))) %>%
    mutate(vars_adj = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AEROD_v",ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "FSDSC",ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"TREFHT","TREFHT"))),
           rep = i
    )
  
  zfi_spatial_df_TREFHT$date <- as_date(zfi_spatial_df_TREFHT$date)
  
  
  preds_plus_T050 <-
    zfi_spatial_df_T050 %>%
    left_join(e3sm_locations, by = c('date', 'loc_id'))  %>%
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
      zfi_contr_wgt = sqrt(weight * (pred_minus_obs_zeroed^2)) - sqrt(weight * (pred_minus_obs^2)),
      rep = i
    ) %>%
    group_by(lat,vars_adj,date,rep) %>%
    summarize(m_zfi = mean(zfi_contr_wgt),s_zfi=sd(zfi_contr_wgt))%>% #taking mean across ESN ensembles
    mutate(model="T050")
  
  
  preds_plus_TREFHT <-
    zfi_spatial_df_TREFHT %>%
    left_join(e3sm_locations, by = c('date', 'loc_id'))  %>%
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
      zfi_contr_wgt = sqrt(weight * (pred_minus_obs_zeroed^2)) - sqrt(weight * (pred_minus_obs^2)),
      rep = i
    )  %>%
    group_by(lat,vars_adj,date,rep) %>%
    summarize(m_zfi = mean(zfi_contr_wgt),s_zfi=sd(zfi_contr_wgt))%>% #taking mean across ESN ensembles
    mutate(model="TREFHT")   
  
  
  preds_plus <-  rbind(preds_plus_T050,preds_plus_TREFHT)
  
  
  
  list(fi_res=fi_res,preds_plus=preds_plus)
}

stopCluster(cl)


fi <- do.call('rbind',lapply(1:length(out),function(x) out[[x]]$fi_res))
preds_plus <- do.call('rbind',lapply(1:length(out),function(x) out[[x]]$preds_plus))

write.csv(fi,file=paste0("../output/gmd_paper/102-e3sm_fullvar.csv"),row.names=FALSE)
write.csv(preds_plus,file=paste0("../output/gmd_paper/102-e3sm_fullvar_preds_plus.csv"),row.names=FALSE)


