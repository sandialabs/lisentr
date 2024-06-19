# This script computes summary statistics for HSW++ plots

library(dplyr)
library(lubridate)
library(tidyr)
library(tidync)
library(listenr)
library(stringr)
library(foreach)
library(doParallel)



# coutnerfactual
fp_cf = paste0("/gpfs/cldera/data/HSW/post/release_011423/counter_factual_latlon/")

strat_h2_file_cf <- list.files(paste0(fp_cf), pattern = 'T050.h2')
surf_h2_file_cf <- list.files(paste0(fp_cf), pattern = 'T1000.h2')

T050_raw_cf = tidync(x = paste0(fp_cf, strat_h2_file_cf))
T1000_raw_cf = tidync(x = paste0(fp_cf, surf_h2_file_cf))

# convert nc to tibble
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


#compute climatologies
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

# compute summary statistics for each ensemble
hsw_list <- list()
out <- for(i in 1:5){
  
  sub_folder1 <- sub_folder_list[[1]][i]
  AOD_h2_file <- list.files(paste0(fp, sub_folder1), pattern = 'AOD.h2')
  strat_h2_file <- list.files(paste0(fp, sub_folder1), pattern = 'T050.h2')
  surf_h2_file <- list.files(paste0(fp, sub_folder1), pattern = 'T1000.h2')
  
  AOD_raw <- tidync(x = paste0(fp, sub_folder1, AOD_h2_file))
  T050_raw = tidync(x = paste0(fp, sub_folder1, strat_h2_file))
  T1000_raw = tidync(x = paste0(fp, sub_folder1, surf_h2_file))
  
# convert nc to tibble
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
  
# compute climatologies
  hsw_clmtg_spatial_means_and_sds2 <-
    hsw_all %>%
    group_by(lon, lat) %>%
    summarise(
      AOD_mean = 0,#mean(AOD, na.rm = TRUE),
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
    filter(time %% 10 == 0) #subsetting data for speed
  

  
  hsw_list[[i]] <- hsw %>%
    summarize(m_strat=mean(T050_stndzd,na.rm=TRUE),m_surf=mean(T1000_stndzd,na.rm=TRUE),m_aod=mean(AOD_stndzd,na.rm=TRUE),.by=c(time)) %>%
    pivot_longer(m_strat:m_aod,names_to='variable',values_to='value') %>%
    mutate(rep=i)
  
}

hsw_out <- do.call('rbind',hsw_list)

write.csv(hsw_out,file=paste0("../output/hsw/108_HSW_means.csv"),row.names=FALSE)

