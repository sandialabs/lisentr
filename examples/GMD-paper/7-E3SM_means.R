# this script computes summary statistics for E3SM

library(dplyr)
library(lubridate)
library(tidyr)
library(tidync)
library(doParallel)
library(foreach)


cl <- makeCluster(5)
registerDoParallel(cl)

# loop over 5 ESM ensembles
out <- foreach(i=1:5,.packages=c('dplyr','lubridate','tidyr','tidync')) %dopar% {
  fp = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,"/post/atm/180x360_aave/ts/monthly/8yr/")
  fp_cf = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,".cf/post/atm/180x360_aave/ts/monthly/8yr/")
  
  
  BURDENSO4_h2_file <- paste0(fp, "BURDENSO4_199101_199812.nc")
  T050_h2_file <- paste0(fp, "T050_199101_199812.nc")
  TREFHT_h2_file <- paste0(fp, "TREFHT_199101_199812.nc")
  AODSO4_h2_file <- paste0(fp, "AODSO4_199101_199812.nc")
  AEROD_v_h2_file <- paste0(fp, "AEROD_v_199101_199812.nc")
  radLW_h2_file <- paste0(fp, "FLNT_199101_199812.nc")
  radSW_h2_file <- paste0(fp, "FSDSC_199101_199812.nc")
  
  
  BURDENSO4_h2_file_cf <- paste0(fp_cf, "BURDENSO4_199101_199812.nc")
  T050_h2_file_cf <- paste0(fp_cf, "T050_199101_199812.nc")
  TREFHT_h2_file_cf <- paste0(fp_cf, "TREFHT_199101_199812.nc")
  AEROD_v_h2_file_cf <- paste0(fp_cf, "AEROD_v_199101_199812.nc")
  AODSO4_h2_file_cf <- paste0(fp_cf, "AODSO4_199101_199812.nc")
  radLW_h2_file_cf <- paste0(fp_cf, "FLNT_199101_199812.nc")
  radSW_h2_file_cf <- paste0(fp_cf, "FSDSC_199101_199812.nc")
  
  #load nc files
  BURDENSO4_raw <- tidync(BURDENSO4_h2_file)
  T050_raw = tidync(T050_h2_file)
  TREFHT_raw = tidync(TREFHT_h2_file)
  AODSO4_raw = tidync(AODSO4_h2_file)
  
  BURDENSO4_raw_cf <- tidync(BURDENSO4_h2_file_cf)
  T050_raw_cf = tidync(T050_h2_file_cf)
  TREFHT_raw_cf = tidync(TREFHT_h2_file_cf)
  AODSO4_raw_cf = tidync(AODSO4_h2_file_cf)
  
  AEROD_v_raw <- tidync(AEROD_v_h2_file)
  FLNT_raw = tidync(radLW_h2_file)
  FSDSC_raw = tidync(radSW_h2_file)
  
  AEROD_v_raw_cf <- tidync(AEROD_v_h2_file_cf)
  FLNT_raw_cf = tidync(radLW_h2_file_cf)
  FSDSC_raw_cf = tidync(radSW_h2_file_cf)
  
  #convert nc files to tibbles
  AEROD_v <- 
    AEROD_v_raw %>% 
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
  
  BURDENSO4 <- 
    BURDENSO4_raw %>% 
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
  
  
  
  AODSO4 <-
    AODSO4_raw %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  # join data
  e3sm_all <- 
    full_join(BURDENSO4, T050,
              by = c("lon", "lat", "time")) %>%
    full_join(TREFHT,
              by = c("lon", "lat", "time")) %>%
    full_join(AODSO4,
              by = c("lon", "lat", "time")) %>%
    full_join(FLNT,
              by = c("lon", "lat", "time")) %>%
    full_join(FSDSC,
              by = c("lon", "lat", "time")) %>%
    full_join(AEROD_v,
              by = c("lon", "lat", "time")) %>%
    dplyr::select(lon, lat, time, everything()) %>%
    mutate(
      date = as_date("1991-01-31") + days(time-min(time)),
      month=month(date)
    )
  
  
  
  #conver nc files to tibbles for counterfactuals
  AEROD_v_cf <- 
    AEROD_v_raw_cf %>% 
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  FLNT_cf <-
    FLNT_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  FSDSC_cf <-
    FSDSC_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon))
  
  
  BURDENSO4_cf <- 
    BURDENSO4_raw_cf %>% 
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
  
  
  AODSO4_cf <-
    AODSO4_raw_cf %>%
    hyper_tibble() %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) 
  
  
  
  
  # joint counterfactuals
  e3sm_all_cf <- 
    full_join(BURDENSO4_cf, T050_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(TREFHT_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(AODSO4_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(FLNT_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(FSDSC_cf,
              by = c("lon", "lat", "time")) %>%
    full_join(AEROD_v_cf,
              by = c("lon", "lat", "time")) %>%
    dplyr::select(lon, lat, time, everything()) %>%
    mutate(
      date = as_date("1991-01-31") + days(time-min(time)),
      month=month(date)
    )
  
  
  
  # compute climatologies
  e3sm_clmtg_spatial_means_and_sds1_cf <-
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
  
  
# compute averages for e3sm variables over time
  timeavg_reg <- e3sm_all %>% 
    select(-c(time,month)) %>%
    pivot_longer(BURDENSO4:AEROD_v,names_to="variable",values_to="value") %>%
    summarize(m_value=mean(value,na.rm=TRUE),s_value=sd(value,na.rm=TRUE),.by=c(date,variable))%>%
    mutate(
      run="regular",
      fake = sin(month(date)*(2*pi/6))
    )
  
  timeavg_cf <- e3sm_all_cf %>% 
    select(-c(time,month)) %>%
    pivot_longer(BURDENSO4:AEROD_v,names_to="variable",values_to="value") %>%
    summarize(m_value=mean(value,na.rm=TRUE),s_value=sd(value,na.rm=TRUE),.by=c(date,variable))%>%
    mutate(
      run="cf",
      fake = sin(month(date)*(2*pi/6))
    )
  
  
  timeavg <- rbind(timeavg_reg,timeavg_cf)
  
  # compute anomalies
  e3sm <-
    left_join(e3sm_all, e3sm_clmtg_spatial_means_and_sds1_cf, 
              by = c("lon", "lat","month")) %>%
    mutate(
      AEROD_v_stndzd = (AEROD_v - AEROD_v_mean) / (AEROD_v_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd),
      FLNT_stndzd = (FLNT - FLNT_mean) / (FLNT_sd),
      FSDSC_stndzd = (FSDSC - FSDSC_mean) / (FSDSC_sd),
      BURDENSO4_stndzd = (BURDENSO4 - BURDENSO4_mean) / (BURDENSO4_sd),
      AODSO4_stndzd = (AODSO4 - AODSO4_mean) / (AODSO4_sd),
      
      date = ceiling_date(date,'month')
    ) 
  
  e3sm2 <-
    left_join(e3sm_all_cf, e3sm_clmtg_spatial_means_and_sds1_cf, 
              by = c("lon", "lat","month")) %>%
    mutate(
      AEROD_v_stndzd = (AEROD_v - AEROD_v_mean) / (AEROD_v_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd),
      FLNT_stndzd = (FLNT - FLNT_mean) / (FLNT_sd),
      FSDSC_stndzd = (FSDSC - FSDSC_mean) / (FSDSC_sd),
      BURDENSO4_stndzd = (BURDENSO4 - BURDENSO4_mean) / (BURDENSO4_sd),
      AODSO4_stndzd = (AODSO4 - AODSO4_mean) / (AODSO4_sd),
      
      date = ceiling_date(date,'month')
    ) 
  
  
  
  
  timeavg_reg2 <- e3sm %>% 
    pivot_longer(AEROD_v_stndzd:AODSO4_stndzd,names_to="variable",values_to="value") %>%
    summarize(m_value=mean(value,na.rm=TRUE),s_value=sd(value,na.rm=TRUE),.by=c(date,variable))%>%
    mutate(
      run="regular",
      fake = sin(month(date)*(2*pi/6))
    )
  
  timeavg_cf2 <- e3sm2 %>% 
    pivot_longer(AEROD_v_stndzd:AODSO4_stndzd,names_to="variable",values_to="value") %>%
    summarize(m_value=mean(value,na.rm=TRUE),s_value=sd(value,na.rm=TRUE),.by=c(date,variable))%>%
    mutate(
      run="cf",
      fake = sin(month(date)*(2*pi/6))
    )
  
  
  timeavg2 <- rbind(timeavg_reg2,timeavg_cf2) %>% mutate(rep=i)
  timeavg2
}


full <- do.call("rbind",out)
write.csv(full,file="../output/e3sm/84-E3SM_SNR.csv",row.names=FALSE)

stopCluster(cl)




