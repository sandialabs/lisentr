
# This script runs a hyperparameter search for ESN on the hsw data




library(dplyr)
library(lubridate)
library(tidyr)
library(tidync)
library(listenr)
library(stringr)
library(foreach)
library(doParallel)

# time lag 
m= 5

# number of ESN ensembles
n_ensembles = 10


# number of month blocks for FI
nblocks = 3




#set up parallelization
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
  
# extract nc files into tibble form
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
  
  # join all data
  hsw_all <- 
    full_join(AOD, T050,
              by = c("lon", "lat", "time")) %>%
    full_join(T1000,
              by = c("lon", "lat", "time")) %>%
    dplyr::select(lon, lat, time, everything()) 
  
  # compute climatologies by lat/lon
  hsw_clmtg_spatial_means_and_sds1 <-
    hsw_all %>%
    group_by(lon, lat) %>%
    summarise(
      AOD_mean = mean(AOD, na.rm = TRUE),
      T050_mean = mean(T050, na.rm = TRUE),
      T1000_mean = mean(T1000, na.rm = TRUE),
      AOD_sd = sd(AOD, na.rm = TRUE),
      T050_sd = sd(T050, na.rm = TRUE),
      T1000_sd = sd(T1000, na.rm = TRUE),
      .groups = "drop"
    )
  
  # calculate anomalies
  hsw <-
    left_join(hsw_all, hsw_clmtg_spatial_means_and_sds1, 
              by = c("lon", "lat")) %>%
    mutate(
      AOD_stndzd = (AOD - AOD_mean) / (AOD_sd),
      T050_stndzd = (T050 - T050_mean) / (T050_sd),
      T1000_stndzd = (T1000 - T1000_mean) / (T1000_sd)
    ) %>%
    filter(time %% 10 == 0) #subsetting data for speed
  
  
  
  hsw_train <- hsw %>%
    mutate(data=ifelse(time > 800, "test","train"))#%>% inner_join(complete_latlon)
  
  
  # get order of lat/lon for cosine weights
  a <- hsw_train %>%
    filter(data=="train") %>%
    select(time, lon, lat, T1000_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = time,
      values_from = T1000_stndzd
    ) 
  
  lat_flat <- a$lat
  weights <- cos(pi/180*lat_flat)
  
  # create training/test matrices for ESN
  train_mat_AOD <-
    hsw_train %>%
    filter(data=="train") %>%
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
    filter(data=="train") %>%
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
    filter(data=="train") %>%
    select(time, lon, lat, T1000_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = time,
      values_from = T1000_stndzd
    ) %>%
    select(-lon,-lat) %>%
    t()
  
  test_mat_AOD <-
    hsw_train %>%
    filter(data=="test") %>%
    select(time, lon, lat, AOD_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = time,
      values_from = AOD_stndzd
    ) %>%
    select(-lon,-lat) %>%
    t()
  
  test_mat_T050 <-
    hsw_train %>%
    filter(data=="test") %>%
    select(time, lon, lat, T050_stndzd) %>%
    pivot_wider(
      id_cols = c(lon, lat),
      names_from = time,
      values_from = T050_stndzd
    ) %>%
    select(-lon,-lat) %>%
    t()
  
  test_mat_T1000 <-
    hsw_train %>%
    filter(data=="test") %>%
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
  x_T050_te = cbind(test_mat_AOD %*% eofs_AOD$phi, test_mat_T050 %*% eofs_T050$phi)
  
  
  # for convenience
  y_train_mat_T050 <- train_mat_T050#
  phi_train_T050 <- eofs_T050$phi 
  
  x_T1000 = cbind(eofs_AOD$train, eofs_T1000$train)
  y_T1000 = eofs_T1000$train
  x_T1000_te = cbind(test_mat_AOD %*% eofs_AOD$phi, test_mat_T1000 %*% eofs_T1000$phi)
  
  # for convenience
  y_train_mat_T1000 <- train_mat_T1000#T1000 #train_mat_T1000 
  phi_train_T1000 <- eofs_T1000$phi #eofs_T050$phi
  
  t_T050 = rownames(x_T050)
  t_T1000 = rownames(x_T1000)
  
  t_T050_te = rownames(x_T050_te)
  t_T1000_te = rownames(x_T1000_te)
  
  # specify hyperparameters
  tuning_params <-
    list(
      nensm = 5,
      cores = 40,
      tau = 1,
      m = 3,
      tau_emb = 1,
      nh = c(25, 50, 100, 200),
      U_width = c(0.1, 0.5),
      W_width = c(0.1, 0.5),
      U_pi = c(0.1, 0.5),
      W_pi = c(0.1, 0.5),
      nu = c(0.1, 0.5),
      reg_par = c(.5,5,50),
      add_quad = c(FALSE),
      internal_scaling = "joint",
      seed = 20231002
    )
  
  # hyperparameter search
  
  # T050
  start <- Sys.time()
  hs_pcs10_T050 <- hyperparameter_search(
    x = x_T050,
    y = y_T050,
    t = as.character(t_T050),
    x_test = x_T050_te,
    t_test = t_T050_te,
    phi = phi_train_T050,
    obs_train = train_mat_T050,
    obs_test = test_mat_T050,
    tau = tuning_params$tau,
    m = tuning_params$m,
    tau_emb = tuning_params$tau_emb,
    nh = tuning_params$nh,
    U_width = tuning_params$U_width,
    W_width = tuning_params$W_width,
    U_pi = tuning_params$U_pi,
    W_pi = tuning_params$W_pi,
    nu = tuning_params$nu,
    reg_par = tuning_params$reg_par,
    add_quad = tuning_params$add_quad,
    internal_scaling = tuning_params$internal_scaling,
    seed = tuning_params$seed,
    n_ensembles = tuning_params$nensm,
    cores = tuning_params$cores,
    weights = weights
  )
  finish <- Sys.time()
  finish - start
  
  #T1000
  start <- Sys.time()
  hs_pcs10_T1000 <- hyperparameter_search(
    x = x_T1000,
    y = y_T1000,
    t = as.character(t_T1000),
    x_test = x_T1000_te,
    t_test = t_T1000_te,
    phi = phi_train_T1000,
    obs_train = train_mat_T1000,
    obs_test = test_mat_T1000,
    tau = tuning_params$tau,
    m = tuning_params$m,
    tau_emb = tuning_params$tau_emb,
    nh = tuning_params$nh,
    U_width = tuning_params$U_width,
    W_width = tuning_params$W_width,
    U_pi = tuning_params$U_pi,
    W_pi = tuning_params$W_pi,
    nu = tuning_params$nu,
    reg_par = tuning_params$reg_par,
    add_quad = tuning_params$add_quad,
    internal_scaling = tuning_params$internal_scaling,
    seed = tuning_params$seed,
    n_ensembles = tuning_params$nensm,
    cores = tuning_params$cores,
    weights = weights
  )
  
    list(hs_pcs10_T050$output,hs_pcs10_T1000$output)
}

stopCluster(cl)

#combine results
t050_te = Reduce(`+`, (lapply(1:5,function(x) out[[x]][[1]]$mRMSE_te))) / length(out)
t1000_te = Reduce(`+`, (lapply(1:5,function(x) out[[x]][[2]]$mRMSE_te))) / length(out)

t050_tr = Reduce(`+`, (lapply(1:5,function(x) out[[x]][[1]]$mRMSE_tr))) / length(out)
t1000_tr = Reduce(`+`, (lapply(1:5,function(x) out[[x]][[2]]$mRMSE_tr))) / length(out)

res_T050 = data.frame(out[[1]][[1]][,1:11],mRMSE_te=t050_te,mRMSE_tr=t050_tr)
res_T1000 = data.frame(out[[1]][[2]][,1:11],mRMSE_te=t1000_te,mRMSE_tr=t1000_tr)

write.csv(res_T050,file=paste0("../output/gmd_paper/105_HSW_hs_strat_20pcs.csv"),row.names=FALSE)
write.csv(res_T1000,file=paste0("../output/gmd_paper/105_HSW_hs_surf_20pcs.csv"),row.names=FALSE)


