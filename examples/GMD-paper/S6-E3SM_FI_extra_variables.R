# This script trains an ESN and FI to the E3SM fullvar 8 year run data with BURDENSO4 and AODSO4 extra



library(dplyr)
library(lubridate)
library(tidyr)
library(tidync)
library(listenr)
library(stringr)
library(foreach)
library(doParallel)


# number of ESN ensembles
n_ensembles = 10

#m parameter
m= 3

# number of month blocks for FI
nblocks = 4
m=3
cl <- makeCluster(5)
registerDoParallel(cl)

#loop over 5 ensembles
out <- foreach(i=1:5,.packages=c('listenr','dplyr','lubridate','tidyr','tidync','stringr')) %dopar% {
  
  fp = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,"/post/atm/180x360_aave/ts/monthly/8yr/")
  fp_cf = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,".cf/post/atm/180x360_aave/ts/monthly/8yr/")
  
  AEROD_v_h2_file <- paste0(fp, "AEROD_v_199101_199812.nc")
  T050_h2_file <- paste0(fp, "T050_199101_199812.nc")
  TREFHT_h2_file <- paste0(fp, "TREFHT_199101_199812.nc")
  FLNT_h2_file <- paste0(fp, "FLNT_199101_199812.nc")
  FSDSC_h2_file <- paste0(fp, "FSDSC_199101_199812.nc")
  BURDENSO4_h2_file <- paste0(fp, "BURDENSO4_199101_199812.nc")
  AODSO4_h2_file <- paste0(fp, "AODSO4_199101_199812.nc")
  
  AEROD_v_raw <- tidync(AEROD_v_h2_file)
  T050_raw = tidync(T050_h2_file)
  TREFHT_raw = tidync(TREFHT_h2_file)
  FLNT_raw = tidync(FLNT_h2_file)
  FSDSC_raw = tidync(FSDSC_h2_file)
  BURDENSO4_raw <- tidync(BURDENSO4_h2_file)
  AODSO4_raw = tidync(AODSO4_h2_file)
  
  #coutnerfactual for baseline
  AEROD_v_h2_file_cf <- paste0(fp_cf, "AEROD_v_199101_199812.nc")
  T050_h2_file_cf <- paste0(fp_cf, "T050_199101_199812.nc")
  TREFHT_h2_file_cf <- paste0(fp_cf, "TREFHT_199101_199812.nc")
  FLNT_h2_file_cf <- paste0(fp_cf, "FLNT_199101_199812.nc")
  FSDSC_h2_file_cf <- paste0(fp_cf, "FSDSC_199101_199812.nc")
  BURDENSO4_h2_file_cf <- paste0(fp_cf, "BURDENSO4_199101_199812.nc")
  AODSO4_h2_file_cf <- paste0(fp_cf, "AODSO4_199101_199812.nc")
  
  AEROD_v_raw_cf <- tidync(AEROD_v_h2_file_cf)
  T050_raw_cf = tidync(T050_h2_file_cf)
  TREFHT_raw_cf = tidync(TREFHT_h2_file_cf)
  FLNT_raw_cf = tidync(FLNT_h2_file_cf)
  FSDSC_raw_cf = tidync(FSDSC_h2_file_cf)
  BURDENSO4_raw_cf <- tidync(BURDENSO4_h2_file_cf)
  AODSO4_raw_cf = tidync(AODSO4_h2_file_cf)
  
  #conver nc files to tibble
  AEROD_v <- 
    AEROD_v_raw %>% 
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) # %>%
  
  T050 <-
    T050_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  
  TREFHT <-
    TREFHT_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  
  FLNT <-
    FLNT_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  
  FSDSC <-
    FSDSC_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) #%>%
  
  BURDENSO4 <- 
    BURDENSO4_raw %>% 
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon))  
  
  
  AODSO4 <-
    AODSO4_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon))
  
  #join data
  e3sm_all <- 
    full_join(AEROD_v, T050,
              by = c("lon", "lat", "time")) %>%
    full_join(TREFHT,
              by = c("lon", "lat", "time")) %>%
    full_join(FLNT,
              by = c("lon", "lat", "time")) %>%
    full_join(FSDSC,
              by = c("lon", "lat", "time")) %>%
    full_join(BURDENSO4,
              by = c("lon", "lat", "time")) %>%
    full_join(AODSO4,
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
  
  BURDENSO4_cf <- 
    BURDENSO4_raw_cf %>% 
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  AODSO4_cf <-
    AODSO4_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon))
  
  e3sm_all_cf <- 
    full_join(AEROD_v_cf, T050_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(TREFHT_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(FLNT_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(FSDSC_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(BURDENSO4_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(AODSO4_cf,
              by = c("lon", "lat", "time")) %>%
    dplyr::select(lon, lat, time, everything()) %>%
    mutate(
      date = as_date("1991-01-31") + days(time-min(time)),
      month=month(date)
    )
  
  
  
  # climatologies
  e3sm_clmtg_spatial_means_and_sds1 <-
    e3sm_all_cf %>%
    group_by(lon, lat,month) %>%
    summarise(
      AEROD_v_mean = mean(AEROD_v, na.rm = TRUE),
      T050_mean = mean(T050, na.rm = TRUE),
      TREFHT_mean = mean(TREFHT, na.rm = TRUE),
      FLNT_mean = mean(FLNT, na.rm = TRUE),
      FSDSC_mean = mean(FSDSC, na.rm = TRUE),
      BURDENSO4_mean = mean(BURDENSO4, na.rm = TRUE),
      AODSO4_mean = mean(AODSO4, na.rm = TRUE),
      
      AEROD_v_sd = sd(AEROD_v, na.rm = TRUE),
      T050_sd = sd(T050, na.rm = TRUE),
      TREFHT_sd = sd(TREFHT, na.rm = TRUE),
      FLNT_sd = sd(FLNT, na.rm = TRUE),
      FSDSC_sd = sd(FSDSC, na.rm = TRUE),
      BURDENSO4_sd = sd(BURDENSO4, na.rm = TRUE),
      AODSO4_sd = sd(AODSO4, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  # anomalies
  e3sm <-
    left_join(e3sm_all, e3sm_clmtg_spatial_means_and_sds1, 
              by = c("lon", "lat","month")) %>%
    mutate(
      AEROD_v_stndzd = (AEROD_v - AEROD_v_mean) / (AEROD_v_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd),
      FLNT_stndzd = (FLNT - FLNT_mean) / (FLNT_sd),
      FSDSC_stndzd = (FSDSC - FSDSC_mean) / (FSDSC_sd),
      BURDENSO4_stndzd = (BURDENSO4 - BURDENSO4_mean) / (BURDENSO4_sd),
      AODSO4_stndzd = (AODSO4 - AODSO4_mean) / (AODSO4_sd)
      
    )
  
  
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
  
  #training matrices
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
  
  train_mat_BURDENSO4 <-
    e3sm_train %>%
    select(date, lon, lat, BURDENSO4_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = BURDENSO4_stndzd
    ) %>%
    select(-lon,-lat) %>%
    t()
  
  train_mat_AODSO4 <-
    e3sm_train %>%
    select(date, lon, lat, AODSO4_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = AODSO4_stndzd
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
  eofs_BURDENSO4 = compute_eofs(Ztrain = train_mat_BURDENSO4, n_eofs = n_eofs)
  eofs_AODSO4 = compute_eofs(Ztrain = train_mat_AODSO4, n_eofs = n_eofs)
  
  
  
  # Model matrices
  x_T050 = cbind(eofs_AEROD_v$train, eofs_FLNT$train, eofs_BURDENSO4$train, eofs_AODSO4$train, eofs_T050$train)
  y_T050 = eofs_T050$train
  
  # so you don't have to change for compute_fi 
  y_train_mat_T050 <- train_mat_T050#
  phi_train_T050 <- eofs_T050$phi 
  
  x_TREFHT = cbind(eofs_AEROD_v$train, eofs_FSDSC$train, eofs_BURDENSO4$train, eofs_AODSO4$train, eofs_TREFHT$train)
  y_TREFHT = eofs_TREFHT$train
  
  y_train_mat_TREFHT <- train_mat_TREFHT 
  phi_train_TREFHT <- eofs_TREFHT$phi 
  
  
  
  t = rownames(x_T050)
  
  #train ESN
  esn_T050 <-
    fit_Eesn(
      x = x_T050,
      y = y_T050,
      t = as.character(t),
      tau = 1,
      m = m,
      tau_emb = 1,
      nh = 50,
      nu = 0.1,
      W_pi = 0.5,
      U_pi = 0.5,
      add_quad = FALSE,
      internal_scaling = "joint",
      seed = 20230413,
      reg_par = 5,
      n_ensembles = n_ensembles
    )
  
  esn_TREFHT <-
    fit_Eesn(
      x = x_TREFHT,
      y = y_TREFHT,
      t = as.character(t),
      tau = 1,
      m = m,
      tau_emb = 1,
      nh = 200,
      nu = 0.1,
      add_quad = FALSE,
      internal_scaling = "joint",
      seed = 20230413,
      reg_par = 5,
      n_ensembles = n_ensembles
    )
  
  
  
  
  
  # Compute ZFI
  zfi_wo_retrain_spatial_T050 <-
    compute_fi(
      model = esn_T050,
      type="zfi",
      var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs),(2*n_eofs + 1):(3 * n_eofs),(3*n_eofs + 1):(4 * n_eofs), (4*n_eofs + 1):(5 * n_eofs)), 
      y_spatial = y_train_mat_T050,
      phi = phi_train_T050,
      seed = 20230411,
      blockSize = nblocks,
      return_adj_preds=TRUE,
      weights=weights
    )
  
  zfi_wo_retrain_spatial_TREFHT <-
    compute_fi(
      model = esn_TREFHT,
      type="zfi",
      var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs),(2*n_eofs + 1):(3 * n_eofs),(3*n_eofs + 1):(4 * n_eofs), (4*n_eofs + 1):(5 * n_eofs)),
      y_spatial = y_train_mat_TREFHT,
      phi = phi_train_TREFHT,
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
      variable = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AEROD_v",
                        ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "FLNT",
                               ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"BURDENSO4",
                                      ifelse(vars_adj==paste((3*n_eofs+1):(4*n_eofs),collapse=','),"AODSO4","T050")))), 
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
      variable = ifelse(vars_adj == paste(1:n_eofs,collapse=","), "AEROD_v",
                        ifelse(vars_adj == paste((n_eofs+1):(2*n_eofs),collapse=','), "FSDSC",
                               ifelse(vars_adj==paste((2*n_eofs+1):(3*n_eofs),collapse=','),"BURDENSO4",
                                      ifelse(vars_adj==paste((3*n_eofs+1):(4*n_eofs),collapse=','),"AODSO4","TREFHT")))), 
      t_forecasted = as_date(t_forecasted),
      t_adj = as_date(t_adj)
    ) %>%
    mutate(month_diff = interval(floor_date(as_date(t_adj),'month'),floor_date(as_date(t_forecasted),'month')) %/% months(1)) %>%
    filter(month_diff == 1)   
  
  
  # Feature Importance Results 
  
  full_zfi_wo_retrain_spatial <- rbind(zfi_res_T050 %>% mutate(model="T050"),zfi_res_TREFHT %>% mutate(model="TREFHT"))
  fi_res <- full_zfi_wo_retrain_spatial %>% 
    mutate(rep=i)
  
  
  
  list(fi_res=fi_res)
}

stopCluster(cl)


fi <- do.call('rbind',lapply(1:length(out),function(x) out[[x]]$fi_res))

write.csv(fi,file=paste0("../output/gmd_paper/104h-e3sm_fullvar.csv"),row.names=FALSE)


