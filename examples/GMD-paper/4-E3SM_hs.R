#  hyperparameter search on E3SM ensemble 1

library(dplyr)
library(lubridate)
library(tidyr)
library(tidync)
library(listenr)
library(stringr)
library(foreach)
library(doParallel)

i=1
fp = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,"/post/atm/180x360_aave/ts/monthly/8yr/")
fp_cf = paste0("/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/fullvar/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens",i,".cf/post/atm/180x360_aave/ts/monthly/8yr/")


AEROD_v_h2_file <- paste0(fp, "AEROD_v_199101_199812.nc")
T050_h2_file <- paste0(fp, "T050_199101_199812.nc")
TREFHT_h2_file <- paste0(fp, "TREFHT_199101_199812.nc")
FLNT_h2_file <- paste0(fp, "FLNT_199101_199812.nc")
FSDSC_h2_file <- paste0(fp, "FSDSC_199101_199812.nc")

# load nc files
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

# load nc files for counterfactual
AEROD_v_raw_cf <- tidync(AEROD_v_h2_file_cf)
T050_raw_cf = tidync(T050_h2_file_cf)
TREFHT_raw_cf = tidync(TREFHT_h2_file_cf)
FLNT_raw_cf = tidync(FLNT_h2_file_cf)
FSDSC_raw_cf = tidync(FSDSC_h2_file_cf)



# convert nc files to tibble
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

# join data
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



#counterfactuals nc to tibbles
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
  left_join(e3sm_all, e3sm_clmtg_spatial_means_and_sds1, 
            by = c("lon", "lat","month")) %>%
  mutate(
    AEROD_v_stndzd = (AEROD_v - AEROD_v_mean) / (AEROD_v_sd),
    T050_stndzd = (T050 - T050_mean) / (T050_sd),
    TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd),
    FLNT_stndzd = (FLNT - FLNT_mean) / (FLNT_sd),
    FSDSC_stndzd = (FSDSC - FSDSC_mean) / (FSDSC_sd)
    
  )%>%
  mutate(data=ifelse(year > 1995,"test","train"))


#-------------------
# remove NAs for AOD

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

# training/test matrices 
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

#test data
e3sm_test <- e3sm %>%
  inner_join(complete_latlon) %>%
  filter(data=="test")


test_mat_AEROD_v <-
  e3sm_test %>%
  select(date, lon, lat, AEROD_v_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = AEROD_v_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

test_mat_T050 <-
  e3sm_test %>%
  select(date, lon, lat, T050_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = T050_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

test_mat_TREFHT <-
  e3sm_test %>%
  select(date, lon, lat, TREFHT_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = TREFHT_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

test_mat_FLNT <-
  e3sm_test %>%
  select(date, lon, lat, FLNT_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = FLNT_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

test_mat_FSDSC <-
  e3sm_test %>%
  select(date, lon, lat, FSDSC_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = FSDSC_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

# EOFS
n_eofs <- 20
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
phi_train_T050 <- eofs_T050$phi 

x_TREFHT = cbind(eofs_AEROD_v$train, eofs_FSDSC$train, eofs_TREFHT$train)
y_TREFHT = eofs_TREFHT$train

x_TREFHT_te = cbind(test_mat_AEROD_v %*% eofs_AEROD_v$phi, test_mat_FSDSC %*% eofs_FSDSC$phi, test_mat_TREFHT %*% eofs_TREFHT$phi)

# for convenience 
phi_train_TREFHT <- eofs_TREFHT$phi 

t_T050 = rownames(x_T050)
t_T050_te = rownames(x_T050_te)

t_TREFHT = rownames(x_TREFHT)
t_TREFHT_te = rownames(x_TREFHT_te)

# hyperparameter list
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
#T050
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

write.csv(hs_pcs10_T050$output, file = "../output/gmd_paper/94e-hs_20pcs_T050_red.csv", row.names = FALSE)


# TREFHT
# 
start <- Sys.time()
hs_pcs10_TREFHT <- hyperparameter_search(
  x = x_TREFHT,
  y = y_TREFHT,
  t = as.character(t_TREFHT),
  x_test = x_TREFHT_te,
  t_test = t_TREFHT_te,
  phi = phi_train_TREFHT,
  obs_train = train_mat_TREFHT,
  obs_test = test_mat_TREFHT,
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

write.csv(hs_pcs10_TREFHT$output, file = "../output/gmd_paper/94e-hs_20pcs_TREFHT_red.csv", row.names = FALSE)

