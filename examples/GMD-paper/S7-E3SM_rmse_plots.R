# This script trains an ESN and computes RMSEs



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
  
  
  #convert nc to tibble
  
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
  
  
  
  # compute climatologies
  e3sm_clmtg_spatial_means_and_sds <-
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
  
  #compute anamolies
  e3sm <-
    left_join(e3sm_all, e3sm_clmtg_spatial_means_and_sds, 
              by = c("lon", "lat","month")) %>%
    mutate(
      AEROD_v_stndzd = (AEROD_v - AEROD_v_mean) / (AEROD_v_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd),
      FLNT_stndzd = (FLNT - FLNT_mean) / (FLNT_sd),
      FSDSC_stndzd = (FSDSC - FSDSC_mean) / (FSDSC_sd)
      
    )
  
  
  # perform repeated cross validation to evaluate predictive cpaability across increasingly large training data sets
  preds_T050 <- list()
  count <- 1
  for(yr in 2:6){
    e3sm = e3sm %>% mutate(data=ifelse(year > 1991 + yr,"test","train"))
    
    b1=e3sm %>% group_by(lon,lat) %>% summarize(n=sum(is.na(AEROD_v)),n2=sum(is.na(FSDSC)),n3=n+n2)
    
    complete_latlon <- b1[b1$n3==0,1:2]
    
    #---------------
    
    
    e3sm_train <- e3sm %>%
      inner_join(complete_latlon)%>%
      filter(data=="train")
    
    
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
    
    train_mat_AEROD_v <-
      e3sm_train %>%
      dplyr::select(date, lon, lat, AEROD_v_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = AEROD_v_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    
    
    train_mat_FLNT <-
      e3sm_train %>%
      dplyr::select(date, lon, lat, FLNT_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = FLNT_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    train_mat_FSDSC <-
      e3sm_train %>%
      dplyr::select(date, lon, lat, FSDSC_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = FSDSC_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    train_mat_T050 <-
      e3sm_train %>%
      dplyr::select(date, lon, lat, T050_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = T050_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    train_mat_TREFHT <-
      e3sm_train %>%
      dplyr::select(date, lon, lat, TREFHT_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = TREFHT_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    #test data
    e3sm_test <- e3sm %>%
      inner_join(complete_latlon) %>%
      filter(data=="test")
    
    
    test_mat_AEROD_v <-
      e3sm_test %>%
      dplyr::select(date, lon, lat, AEROD_v_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = AEROD_v_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    test_mat_T050 <-
      e3sm_test %>%
      dplyr::select(date, lon, lat, T050_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = T050_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    test_mat_TREFHT <-
      e3sm_test %>%
      dplyr::select(date, lon, lat, TREFHT_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = TREFHT_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    test_mat_FLNT <-
      e3sm_test %>%
      dplyr::select(date, lon, lat, FLNT_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = FLNT_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    test_mat_FSDSC <-
      e3sm_test %>%
      dplyr::select(date, lon, lat, FSDSC_stndzd) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = FSDSC_stndzd
      ) %>%
      dplyr::select(-lon,-lat) %>%
      t()
    
    
    ## EOFs
    
    n_eofs <- 20
    n_eofs_strat <- 20
    n_eofs_surf = 20
    
    eofs_AEROD_v = compute_eofs(Ztrain = train_mat_AEROD_v, n_eofs = n_eofs)
    eofs_FSDSC = compute_eofs(Ztrain = train_mat_FSDSC, n_eofs = n_eofs)
    eofs_FLNT = compute_eofs(Ztrain = train_mat_FLNT, n_eofs = n_eofs)
    eofs_T050 = compute_eofs(Ztrain = train_mat_T050, n_eofs = n_eofs)
    eofs_TREFHT = compute_eofs(Ztrain = train_mat_TREFHT, n_eofs = n_eofs)
    
    
    
    # Model matrices
    x_T050 = cbind(eofs_AEROD_v$train, eofs_FLNT$train, eofs_T050$train)
    y_T050 = eofs_T050$train
    
    x_T050_te = cbind(test_mat_AEROD_v %*% eofs_AEROD_v$phi, test_mat_FLNT %*% eofs_FLNT$phi, test_mat_T050 %*% eofs_T050$phi)
    
    # for convenience
    y_train_mat_T050 <- train_mat_T050#
    phi_train_T050 <- eofs_T050$phi
    
    x_TREFHT = cbind(eofs_AEROD_v$train[,1:n_eofs_surf], eofs_FSDSC$train[,1:n_eofs_surf], eofs_TREFHT$train[,1:n_eofs_surf])
    y_TREFHT = eofs_TREFHT$train[,1:n_eofs_surf]
    x_TREFHT_te = cbind(test_mat_AEROD_v %*% eofs_AEROD_v$phi[,1:n_eofs_surf], 
                        test_mat_FSDSC %*% eofs_FSDSC$phi[,1:n_eofs_surf], test_mat_TREFHT %*% eofs_TREFHT$phi[,1:n_eofs_surf])
    
    y_train_mat_TREFHT <- train_mat_TREFHT
    phi_train_TREFHT <- eofs_TREFHT$phi[,1:n_eofs_surf] 
    
    
    
    t = rownames(x_T050)
    t_te = rownames(x_T050_te)
    
    #train ESN
    esn_T050 <-
      fit_Eesn(
        x = x_T050,
        y = y_T050,
        t = as.character(t),
        tau = 1,
        m = m,
        tau_emb = 1,
        nh = 200,
        nu = 0.1,
        W_pi = 0.5,
        U_pi = 0.5,
        add_quad = FALSE,
        internal_scaling = "joint",
        seed = 20230413,
        reg_par = 50,
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
    
    
    
    e3sm_locations <- 
      e3sm_train %>%
      dplyr::select(date, lon, lat, T050) %>%
      pivot_wider(
        id_cols = c(lon, lat),
        names_from = date,
        values_from = T050
      ) %>%
      mutate(loc_id = 1:n()) %>%
      dplyr::select(loc_id, lon, lat) #%>%
    #left_join(e3sm_train,by=c('lat','lon'))
    
    preds_list_T050 <- list()
    preds_list_TREFHT <- list()
    
    # loop over esn ensembles to get predictions to then compute RMSEs using EESN
    for(j in 1:n_ensembles){
      #get training data on spatial scale
      preds_tr_T050 <- 
        data.frame(predict_esn(esn_T050[[j]])$preds_ins %*% t(phi_train_T050)) %>%
        tibble::rownames_to_column(var = "t_forecasted") %>%
        pivot_longer(
          cols = -c(t_forecasted),
          names_to = "loc_id",
          values_to = "T050_stndzd_pred"
        ) %>%
        mutate(loc_id = as.numeric(str_remove(loc_id, "X"))) %>%
        mutate(t_forecasted = as_date(t_forecasted),date=t_forecasted,month=month(date)) %>%
        left_join(e3sm_locations, by = "loc_id") %>%
        left_join(e3sm_clmtg_spatial_means_and_sds,by=c("lat","lon","month")) %>%
        mutate(type="Train") %>%
        left_join(e3sm %>% dplyr::select(date,lat,lon,T050_stndzd),by=c("date","lat","lon")) 
      
      #get testing data on spatial scale
      preds_te_T050 <- 
        data.frame(predict_esn(esn_T050[[j]],x_new = x_T050_te,t_new=t_te)$preds_new %*% t(phi_train_T050)) %>%
        tibble::rownames_to_column(var = "t_forecasted") %>%
        pivot_longer(
          cols = -c(t_forecasted),
          names_to = "loc_id",
          values_to = "T050_stndzd_pred"
        ) %>%
        mutate(loc_id = as.numeric(str_remove(loc_id, "X"))) %>%
        mutate(t_forecasted = as_date(substr(t_forecasted,1,10)),date=t_forecasted,month=month(date)) %>%
        left_join(e3sm_locations, by = "loc_id") %>%
        left_join(e3sm_clmtg_spatial_means_and_sds,by=c("lat","lon","month")) %>%
        mutate(type="Test") %>%
        left_join(e3sm  %>% dplyr::select(date,lat,lon,T050_stndzd),by=c("date","lat","lon")) 
      
      preds_T050_combined <- rbind(preds_tr_T050,preds_te_T050) %>%
        dplyr::select(date,lat,lon,T050_stndzd,T050_stndzd_pred,type) #%>%
      #mutate(year=1991+yr,ensemble=j)
      
      preds_tr_TREFHT <- 
        #data.frame(tr$preds_ins) %*% phi_train) %>%
        data.frame(predict_esn(esn_TREFHT[[j]])$preds_ins %*% t(phi_train_TREFHT)) %>%
        tibble::rownames_to_column(var = "t_forecasted") %>%
        pivot_longer(
          cols = -c(t_forecasted),
          names_to = "loc_id",
          values_to = "TREFHT_stndzd_pred"
        ) %>%
        mutate(loc_id = as.numeric(str_remove(loc_id, "X"))) %>%
        mutate(t_forecasted = as_date(t_forecasted),date=t_forecasted,month=month(date)) %>%
        left_join(e3sm_locations, by = "loc_id") %>%
        left_join(e3sm_clmtg_spatial_means_and_sds,by=c("lat","lon","month")) %>%
        mutate(type="Train") %>%
        left_join(e3sm %>% dplyr::select(date,lat,lon,TREFHT_stndzd),by=c("date","lat","lon")) 
      
      #get testing data on spatial scale
      preds_te_TREFHT <- 
        data.frame(predict_esn(esn_TREFHT[[j]],x_new = x_TREFHT_te,t_new=t_te)$preds_new %*% t(phi_train_TREFHT)) %>%
        tibble::rownames_to_column(var = "t_forecasted") %>%
        pivot_longer(
          cols = -c(t_forecasted),
          names_to = "loc_id",
          values_to = "TREFHT_stndzd_pred"
        ) %>%
        mutate(loc_id = as.numeric(str_remove(loc_id, "X"))) %>%
        mutate(t_forecasted = as_date(substr(t_forecasted,1,10)),date=t_forecasted,month=month(date)) %>%
        left_join(e3sm_locations, by = "loc_id") %>%
        left_join(e3sm_clmtg_spatial_means_and_sds,by=c("lat","lon","month")) %>%
        mutate(type="Test") %>%
        left_join(e3sm  %>% dplyr::select(date,lat,lon,TREFHT_stndzd),by=c("date","lat","lon")) 
      
      preds_TREFHT_combined <- rbind(preds_tr_TREFHT,preds_te_TREFHT) %>%
        dplyr::select(date,lat,lon,TREFHT_stndzd,TREFHT_stndzd_pred,type) #%>%
      #mutate(year=1991+yr,ensemble=j)
      
      
      preds_list_T050[[j]] <- preds_T050_combined
      preds_list_TREFHT[[j]] <- preds_TREFHT_combined
      
      
    }
    
    temp_preds= do.call('rbind',preds_list_T050) %>% 
      summarize(T050_stndzd=mean(T050_stndzd,na.rm=TRUE),T050_pred_stndzd=mean(T050_stndzd_pred,na.rm=TRUE),.by=c(date,lat,lon)) %>%
      mutate(weights=cos(pi/180*lat)) %>% 
      summarize(wrmse=sqrt(mean(weights*(T050_stndzd-T050_pred_stndzd)^2)),.by=c(date)) %>%
      mutate(rep=i,year=1991+yr,model="T050")
    
    temp_preds2= do.call('rbind',preds_list_TREFHT) %>% 
      summarize(TREFHT_stndzd=mean(TREFHT_stndzd,na.rm=TRUE),TREFHT_pred_stndzd=mean(TREFHT_stndzd_pred,na.rm=TRUE),.by=c(date,lat,lon)) %>%
      mutate(weights=cos(pi/180*lat)) %>% 
      summarize(wrmse=sqrt(mean(weights*(TREFHT_stndzd-TREFHT_pred_stndzd)^2)),.by=c(date)) %>%
      mutate(rep=i,year=1991+yr,model="TREFHT")
    
    preds_T050[[count]] = rbind(temp_preds  ,temp_preds2)  
    count = count + 1
  }
  
  
  do.call('rbind',preds_T050)
  
}

stopCluster(cl)


out_df <- do.call('rbind',out)

write.csv(out_df,file=paste0("../output/gmd_paper/103-e3sm_fullvar_rmse.csv"),row.names=FALSE)



