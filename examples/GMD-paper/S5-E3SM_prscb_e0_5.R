#this script trains ESN and computes FI for E3SM prescribed ensembles 0.5x pinatubo

library(dplyr)
library(lubridate)
library(tidyr)
library(tidync)
library(listenr)
library(stringr)
library(doParallel)
library(foreach)



# user selected 
rep_num <- c(1,3,4,6)#c(1:4,6)
size <- "e0.5_cori"
n_ensembles <- 20
nblocks = 4
n_eofs <- 10
m = 3
cl <- makeCluster(5)
registerDoParallel(cl)

#loop over ensembles
out <- foreach(i=1:length(rep_num),.packages=c('listenr','dplyr','lubridate','tidyr','tidync','stringr')) %dopar% {
  rep = rep_num[i]
  fp = paste0("/projects/cldera/data/E3SM/e3sm-cldera/le/v2.LR.WCYCL20TR_01",rep,"1_b1991_",size,"/post/atm/180x360_aave/ts/monthly/5yr/")
  fp_cf = paste0("/projects/cldera/data/E3SM/e3sm-cldera/le/v2.LR.WCYCL20TR_01",rep,"1_b1991_e0.0_cori/post/atm/180x360_aave/ts/monthly/5yr/")
  
  AEROD_v_h2_file <- paste0(fp, "AEROD_v_199101_199512.nc")
  T050_h2_file <- paste0(fp, "T050_199101_199512.nc")
  TREFHT_h2_file <- paste0(fp, "TREFHT_199101_199512.nc")
  radLW_h2_file <- paste0(fp, "FLNT_199101_199512.nc")
  radSW_h2_file <- paste0(fp, "FSNT_199101_199512.nc")
  
  AEROD_v_raw <- tidync(AEROD_v_h2_file)
  T050_raw = tidync(T050_h2_file)
  TREFHT_raw = tidync(TREFHT_h2_file)
  FLNT_raw = tidync(radLW_h2_file)
  FSNT_raw = tidync(radSW_h2_file)
  
  #coutnerfactual for baseline
  AEROD_v_h2_file_cf <- paste0(fp_cf, "AEROD_v_199101_199512.nc")
  T050_h2_file_cf <- paste0(fp_cf, "T050_199101_199512.nc")
  TREFHT_h2_file_cf <- paste0(fp_cf, "TREFHT_199101_199512.nc")
  radLW_h2_file_cf <- paste0(fp_cf, "FLNT_199101_199512.nc")
  radSW_h2_file_cf <- paste0(fp, "FSNT_199101_199512.nc")#paste0(fp_cf, "FSNT_199101_199512.nc")
  
  
  AEROD_v_raw_cf <- tidync(AEROD_v_h2_file_cf)
  T050_raw_cf = tidync(T050_h2_file_cf)
  TREFHT_raw_cf = tidync(TREFHT_h2_file_cf)
  FLNT_raw_cf = tidync(radLW_h2_file_cf)
  FSNT_raw_cf = tidync(radSW_h2_file_cf)
  
  
  #conver to tibble
  
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
  
  FSNT <-
    FSNT_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  # join data
  e3sm_all <- 
    full_join(AEROD_v, T050,
              by = c("lon", "lat", "time")) %>%
    full_join(TREFHT,
              by = c("lon", "lat", "time")) %>%
    full_join(FLNT,
              by = c("lon", "lat", "time")) %>%
    full_join(FSNT,
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
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  T050_cf <-
    T050_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  TREFHT_cf <-
    TREFHT_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  FLNT_cf <-
    FLNT_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  FSNT_cf <-
    FSNT_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  # join cf data
  e3sm_all_cf <- 
    full_join(AEROD_v_cf, T050_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(TREFHT_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(FLNT_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(FSNT_cf,
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
      FSNT_mean = mean(FSNT, na.rm = TRUE),
      AEROD_v_sd = sd(AEROD_v, na.rm = TRUE),
      T050_sd = sd(T050, na.rm = TRUE),
      TREFHT_sd = sd(TREFHT, na.rm = TRUE),
      FLNT_sd = sd(FLNT, na.rm = TRUE),
      FSNT_sd = sd(FSNT, na.rm = TRUE),
      .groups = "drop"
    )
  
  # copmute anomalies
  e3sm <-
    left_join(e3sm_all, e3sm_clmtg_spatial_means_and_sds1, 
              by = c("lon", "lat","month")) %>%
    mutate(
      AEROD_v_stndzd = (AEROD_v - AEROD_v_mean) / (AEROD_v_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd),
      FLNT_stndzd = (FLNT - FLNT_mean) / (FLNT_sd),
      FSNT_stndzd = (FSNT - FSNT_mean) / (FSNT_sd)
      
    )
  
  # remove AOD NAs
  b1=e3sm %>% group_by(lon,lat) %>% summarize(n=sum(is.na(AEROD_v)),n2=sum(is.na(FSNT)),n3=n+n2)
  
  complete_latlon <- b1[b1$n3==0,1:2]
  
  
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
  
  # training matrices
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
  
  train_mat_FSNT <-
    e3sm_train %>%
    select(date, lon, lat, FSNT_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = FSNT_stndzd
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
  
  eofs_AEROD_v = compute_eofs(Ztrain = train_mat_AEROD_v, n_eofs = n_eofs)
  eofs_FSNT = compute_eofs(Ztrain = train_mat_FSNT, n_eofs = n_eofs)
  eofs_FLNT = compute_eofs(Ztrain = train_mat_FLNT, n_eofs = n_eofs)
  eofs_T050 = compute_eofs(Ztrain = train_mat_T050, n_eofs = n_eofs)
  eofs_TREFHT = compute_eofs(Ztrain = train_mat_TREFHT, n_eofs = n_eofs)
  
  
  
  # Model
  
  
  x_T050 = cbind(eofs_AEROD_v$train, eofs_FLNT$train, eofs_T050$train)
  y_T050 = eofs_T050$train
  
  # so you don't have to change for compute_fi
  y_train_mat_T050 <- train_mat_T050 
  phi_train_T050 <- eofs_T050$phi 
  
  x_TREFHT = cbind(eofs_AEROD_v$train, eofs_FSNT$train, eofs_TREFHT$train)
  y_TREFHT = eofs_TREFHT$train
  
  y_train_mat_TREFHT <- train_mat_TREFHT
  phi_train_TREFHT <- eofs_TREFHT$phi 
  
  
  
  t = rownames(x_T050)
  
  
  # train ESN
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
      var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs),(2*n_eofs + 1):(3 * n_eofs)),#,(3*n_eofs + 1):(4 * n_eofs)), 
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
      var_groups = list(1:n_eofs, (n_eofs + 1):(2 * n_eofs),(2*n_eofs + 1):(3 * n_eofs)),#,(3*n_eofs + 1):(4 * n_eofs)), 
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
  
  
  
  
  list(fi_res=fi_res)
}

stopCluster(cl)


fi <- do.call('rbind',lapply(1:length(out),function(x) out[[x]]$fi_res))

write.csv(fi,file=paste0("../output/gmd_paper/104e-e3sm_pscrb_",size,".csv"),row.names=FALSE)

