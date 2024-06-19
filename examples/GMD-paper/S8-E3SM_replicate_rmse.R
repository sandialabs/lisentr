## This script copmutes the RMSE using other climate model replicates as predicitons to use as a baseline to compare against




library(dplyr)
library(lubridate)
library(tidyr)
library(tidync)
library(listenr)
library(stringr)
library(foreach)
library(doParallel)

# m paraeter 
m= 3

# number of ESN ensembles
n_ensembles = 10


# number of month blocks for FI
nblocks = 4

#load in all 5 replicates
e3sm_train <- list()
for(i in 1:5){
  fp = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,"/post/atm/180x360_aave/ts/monthly/8yr/")
  fp_cf = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,".cf/post/atm/180x360_aave/ts/monthly/8yr/")

  
  AEROD_v_h2_file <- paste0(fp, "AEROD_v_199101_199812.nc")
  T050_h2_file <- paste0(fp, "T050_199101_199812.nc")
  TREFHT_h2_file <- paste0(fp, "TREFHT_199101_199812.nc")
  FLNT_h2_file <- paste0(fp, "FLNT_199101_199812.nc")
  FSDSC_h2_file <- paste0(fp, "FSDSC_199101_199812.nc")
  
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
  
  
  AEROD_v_raw_cf <- tidync(AEROD_v_h2_file_cf)
  T050_raw_cf = tidync(T050_h2_file_cf)
  TREFHT_raw_cf = tidync(TREFHT_h2_file_cf)
  FLNT_raw_cf = tidync(FLNT_h2_file_cf)
  FSDSC_raw_cf = tidync(FSDSC_h2_file_cf)
  
  
# convert nc to tibbles
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
      month=month(date),
      year=year(date)
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
  
  
  
  # copmute climatologies
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
  
  # compute anamolies
  e3sm <-
    left_join(e3sm_all, e3sm_clmtg_spatial_means_and_sds1, 
              by = c("lon", "lat","month")) %>%
    mutate(
      AEROD_v_stndzd = (AEROD_v - AEROD_v_mean) / (AEROD_v_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd),
      FLNT_stndzd = (FLNT - FLNT_mean) / (FLNT_sd),
      FSDSC_stndzd = (FSDSC - FSDSC_mean) / (FSDSC_sd)
      
    )%>%
    mutate(data=ifelse(year > 1999,"test","train"),replicate=i)
  
  #remove AOD NAs
  
  b1=e3sm %>% group_by(lon,lat) %>% summarize(n=sum(is.na(AEROD_v)),n2=sum(is.na(FSDSC)),n3=n+n2)
  
  complete_latlon <- b1[b1$n3==0,1:2]
  
  #---------------
  
  
  e3sm_train[[i]] <- e3sm %>%
    inner_join(complete_latlon)%>%
    filter(data=="train")
}




#test data
#test data
test_mat_T050 <- list()
test_mat_TREFHT <- list()
for(i in (1:5)){
  test_mat_T050[[i]] <-
    e3sm_train[[i]] %>%
    dplyr::select(date, lon, lat, T050_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = T050_stndzd
    ) %>%
    dplyr::select(-lon,-lat) %>%
    t()
  
  test_mat_TREFHT[[i]] <-
    e3sm_train[[i]] %>%
    dplyr::select(date, lon, lat, TREFHT_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = TREFHT_stndzd
    ) %>%
    dplyr::select(-lon,-lat) %>%
    t()
  
}


cl <- makeCluster(5)
registerDoParallel(cl)



# compute RMSE for each replicate
out <- foreach(i=1:5,.packages=c('listenr','dplyr','lubridate','tidyr','tidync','stringr')) %dopar% {
  
  # define helper functions
  train_rmse <- function(esn_preds_spatial_T050_tr){
    esn_preds_spatial_df_T050 <-
      data.frame(esn_preds_spatial_T050_tr$preds_ins) %>%
      tibble::rownames_to_column(var = "date") %>%
      pivot_longer(
        cols = -date,
        names_to = "loc_id",
        values_to = "pred"
      ) %>%
      mutate(loc_id = as.integer(str_remove(loc_id, "X")),
             date=as_date(date)) 
    
    preds_plus_T050 <- e3sm_locations  %>%
      inner_join(esn_preds_spatial_df_T050, by = c('date', 'loc_id')) %>%
      mutate(weight = sqrt(cos(pi/180*lat)))    %>%
      group_by(date) %>%
      mutate(total_weight = sum(weight, na.rm = T)) %>%
      ungroup() %>%
      mutate(
        pred_minus_obs = pred - T050_stndzd,
      ) %>%
      mutate(
        pred_minus_obs_wgt = weight * pred_minus_obs,
      )
    
    out = preds_plus_T050 %>%
      summarize(wrmse=sqrt(mean(pred_minus_obs_wgt^2,na.rm=TRUE)),.by=c(date)) %>%
      pull(wrmse)
    
    return(out)
  }
  
  test_rmse <- function(obs,pred){
    obs_df <-
      data.frame(obs) %>%
      tibble::rownames_to_column(var = "date") %>%
      pivot_longer(
        cols = -date,
        names_to = "loc_id",
        values_to = "obs"
      ) %>%
      mutate(loc_id = as.integer(str_remove(loc_id, "X")),
             date=as_date(gsub("\\+ 1", "",date))) 
    
    pred_df <-
      data.frame(pred) %>%
      tibble::rownames_to_column(var = "date") %>%
      pivot_longer(
        cols = -date,
        names_to = "loc_id",
        values_to = "pred"
      ) %>%
      mutate(loc_id = as.integer(str_remove(loc_id, "X")),
             date=as_date(gsub("\\+ 1", "",date))) 
    
    preds_plus_T050 <- e3sm_locations  %>%
      inner_join(obs_df, by = c('date', 'loc_id')) %>%
      inner_join(pred_df,by=c('date','loc_id')) %>% 
      mutate(weight = sqrt(cos(pi/180*lat)))    %>%
      group_by(date) %>%
      mutate(total_weight = sum(weight, na.rm = T)) %>%
      ungroup() %>%
      mutate(
        pred_minus_obs = pred - obs,
      ) %>%
      mutate(
        pred_minus_obs_wgt = weight * pred_minus_obs,
      )
    
    out = preds_plus_T050 %>%
      summarize(wrmse=sqrt(mean(pred_minus_obs_wgt^2,na.rm=TRUE)),.by=c(date)) %>%
      pull(wrmse)
    
    return(out)
  }

  
  
  
  
  # get order of lat/lon for weights
  a <- e3sm_train[[i]] %>%
    dplyr::select(date, lon, lat, TREFHT_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = TREFHT_stndzd
    ) 
  
  lat_flat <- a$lat
  weights <- cos(pi/180*lat_flat)
  
  # create training matrices
  train_mat_AEROD_v <-
    e3sm_train[[i]] %>%
    dplyr::select(date, lon, lat, AEROD_v_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = AEROD_v_stndzd
    ) %>%
    dplyr::select(-lon,-lat) %>%
    t()
  
  
  
  train_mat_FLNT <-
    e3sm_train[[i]] %>%
    dplyr::select(date, lon, lat, FLNT_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = FLNT_stndzd
    ) %>%
    dplyr::select(-lon,-lat) %>%
    t()
  
  train_mat_FSDSC <-
    e3sm_train[[i]] %>%
    dplyr::select(date, lon, lat, FSDSC_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = FSDSC_stndzd
    ) %>%
    dplyr::select(-lon,-lat) %>%
    t()
  
  train_mat_T050 <-
    e3sm_train[[i]] %>%
    dplyr::select(date, lon, lat, T050_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = T050_stndzd
    ) %>%
    dplyr::select(-lon,-lat) %>%
    t()
  
  train_mat_TREFHT <-
    e3sm_train[[i]] %>%
    dplyr::select(date, lon, lat, TREFHT_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = TREFHT_stndzd
    ) %>%
    dplyr::select(-lon,-lat) %>%
    t()
  
  ## EOFs
  
  n_eofs <- 10
  eofs_AEROD_v = compute_eofs(Ztrain = train_mat_AEROD_v, n_eofs = n_eofs)
  eofs_FSDSC = compute_eofs(Ztrain = train_mat_FSDSC, n_eofs = n_eofs)
  eofs_FLNT = compute_eofs(Ztrain = train_mat_FLNT, n_eofs = n_eofs)
  eofs_T050 = compute_eofs(Ztrain = train_mat_T050, n_eofs = n_eofs)
  eofs_TREFHT = compute_eofs(Ztrain = train_mat_TREFHT, n_eofs = n_eofs)
  
  
  
  # Model matrices
  x_T050 = cbind(eofs_AEROD_v$train, eofs_FLNT$train, eofs_T050$train)
  y_T050 = eofs_T050$train
  
  # for convenience
  phi_train_T050 <- eofs_T050$phi 
  
  
  x_TREFHT = cbind(eofs_AEROD_v$train, eofs_FSDSC$train, eofs_TREFHT$train)
  y_TREFHT = eofs_TREFHT$train
  
  phi_train_TREFHT <- eofs_TREFHT$phi 
  
  
  t = rownames(x_T050)
  
  #train ESNs
  esn_T050 <-
    fit_Eesn(
      x = x_T050,
      y = y_T050,
      t = as.character(t),
      tau = 1,
      m = m,
      tau_emb = 1,
      nh = 50,
      add_quad = TRUE,
      internal_scaling = "joint",
      seed = 20230413,
      reg_par = .1,
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
      nh = 50,
      add_quad = TRUE,
      internal_scaling = "joint",
      seed = 20230413,
      reg_par = .1,
      n_ensembles = n_ensembles
    )
  
  e3sm_locations <- 
    e3sm %>%
    dplyr::select(date, lon, lat, T050) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = date,
      values_from = T050
    ) %>%
    mutate(loc_id = 1:n()) %>%
    dplyr::select(loc_id, lon, lat) %>%
    left_join(e3sm,by=c('lat','lon'))
  
  # loop over ESN ensembles to get EESN
  wrmse_T050_list <- list()
  wrmse_TREFHT_list <- list()
  for(j in 1:n_ensembles){
    #Compute ESN predictions on spatial scale:
    esn_preds_spatial_T050_tr <- predict_esn(model = esn_T050[[j]], phi = phi_train_T050)
    esn_preds_spatial_TREFHT_tr <- predict_esn(model = esn_TREFHT[[j]], phi = phi_train_TREFHT)
    
    
    wrmse_T050_list[[j]] <- train_rmse(esn_preds_spatial_T050_tr)
    wrmse_TREFHT_list[[j]] <- train_rmse(esn_preds_spatial_TREFHT_tr)
    
  }
  
  # loop over E3SM ensembles to get out of sample esitmate of RMSE
  wrmse_T050_te_list <- list()
  wrmse_TREFHT_te_list <- list()
  for(k in (1:5)[-i]){
    # need to excldue the first m+1 predictions sisnce we don't get those with model prediction
    wrmse_T050_te_list[[k]] <- test_rmse(train_mat_T050,test_mat_T050[[k]])[-c(1:(m+1))] 
    wrmse_TREFHT_te_list[[k]] <- test_rmse(train_mat_TREFHT,test_mat_TREFHT[[k]])[-c(1:(m+1))]
  }
  
  dates_tr = row.names(train_mat_AEROD_v)[-c(1:(m+1))]
  
  
  wrmse_T050 = rowMeans(do.call('cbind',wrmse_T050_list))
  wrmse_T050_te = rowMeans(do.call('cbind',wrmse_T050_te_list))
  
  wrmse_TREFHT = rowMeans(do.call('cbind',wrmse_TREFHT_list))
  wrmse_TREFHT_te = rowMeans(do.call('cbind',wrmse_TREFHT_te_list))
  
  T050_out = data.frame(train=wrmse_T050,test=wrmse_T050_te) %>% 
    mutate(date=dates_tr,rep=i,model='T050') 
  TREFHT_out = data.frame(train=wrmse_TREFHT,test=wrmse_TREFHT_te) %>% 
    mutate(date=dates_tr,rep=i,model='TREFHT')  
  
  rbind(T050_out,TREFHT_out)
  
}

stopCluster(cl)


out_df <- do.call('rbind',out)

write.csv(out_df,file=paste0("../output/e3sm/95-e3sm_replicate_rmse.csv"),row.names=FALSE)



