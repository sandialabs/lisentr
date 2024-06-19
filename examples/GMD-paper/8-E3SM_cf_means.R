# this script computes summary statistics for e3sm counterfactual

library(dplyr)
library(lubridate)
library(tidyr)
library(tidync)
library(stringr)
library(foreach)
library(doParallel)



cl <- makeCluster(5)
registerDoParallel(cl)

#loop over 5 cf ensembles
out <- foreach(i=1:5,.packages=c('dplyr','lubridate','tidyr','tidync','stringr')) %dopar% {
    fp_cf = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,".cf/post/atm/180x360_aave/ts/monthly/8yr/")
  
  
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
  
  
  #convert nc file to tibble
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
  
  # join data
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
    left_join(e3sm_all_cf, e3sm_clmtg_spatial_means_and_sds1, 
              by = c("lon", "lat","month")) %>%
    mutate(
      AEROD_v_stndzd = (AEROD_v - AEROD_v_mean) / (AEROD_v_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd),
      FLNT_stndzd = (FLNT - FLNT_mean) / (FLNT_sd),
      FSDSC_stndzd = (FSDSC - FSDSC_mean) / (FSDSC_sd)
      
    )
  
  e3sm_out <- e3sm %>%
    group_by(date) %>%
    summarize(m_AEROD_v_stndzd=mean(AEROD_v_stndzd,na.rm=TRUE),
              m_T050_stndzd=mean(T050_stndzd,na.rm=TRUE),
              m_TREFHT_stndzd=mean(TREFHT_stndzd,na.rm=TRUE),
              m_FLNT_stndzd=mean(FLNT_stndzd,na.rm=TRUE),
              m_FSDSC_stndzd=mean(FSDSC_stndzd,na.rm=TRUE)
              ) %>%
    pivot_longer(m_AEROD_v_stndzd:m_FSDSC_stndzd,names_to='variable',values_to='value') %>%
    mutate(rep=i)
    
  e3sm_out
}


outdf <- do.call('rbind',out)
write.csv(outdf,file=paste0("../output/e3sm/74p_e3sm_cf_means.csv"),row.names=FALSE)
