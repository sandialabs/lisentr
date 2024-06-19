# this script performs a hyperparameter search for merra-2

library(dplyr)
library(doParallel)
library(foreach)
library(ggplot2)
library(listenr)
library(lubridate)
library(purrr)
library(tidyr)
library(stringr)


nretrains = 10
n_ensembles = 10

fp = "/projects/cldera/obs_thrust/MERRA-2/"

#load data
TOTEXTTAU_raw = read.csv(paste0(fp, "tavgM_2d_aer_Nx/TOTEXTTAUall_360x180.csv"))
T050_raw = read.csv(paste0(fp, "instM_3d_asm_Np/MERRA2_3dasm_temperature_50mb_360x180.csv"))
TREFHT_raw = read.csv(paste0(fp, "tavg1_2d_slv_Nx_daily/MERRA2_tavg1_2d_slv_monthly_Nx_360x180.csv"))
rad_raw <- read.csv(paste0(fp,"tavgM_2d_rad_Nx/Rad_360x180.csv"))



## Cleaning
TOTEXTTAU <- 
  TOTEXTTAU_raw %>% 
  mutate(date = as_date(date))

rad <- 
  rad_raw %>% 
  mutate(date = as_date(date))

T050 <-
  T050_raw %>% 
  rename("T050" = `T`) %>%
  mutate(date = as_date(date)) %>%
  mutate(day = day(date))

TREFHT <-
  TREFHT_raw %>% 
  rename("TREFHT" = `T`) %>%
  mutate(date = as_date(date)) %>%
  mutate(day = day(date),
         month = month(date),
         year = year(date))



# join data
merra2 <- 
  full_join(TOTEXTTAU, T050 %>% select(-c(month,year,day)), by = c("lon", "lat", "date")) %>%
  full_join(rad %>% filter(date <= as_date("1998-12-31")),by= c("lon", "lat", "date")) %>%
  full_join(TREFHT %>% select(-c(month,year,day)),by= c("lon", "lat", "date")) %>%
  mutate(year=year(date),month=month(date),data=ifelse(year >1995,'test','train')) %>%
  filter(year >= 1991)




# compute climatologies
merra2_clmtg_spatial_means_and_sds <-
  merra2 %>%
  filter(data == "train") %>%
  group_by(month, lon, lat) %>%
  summarise(
    TOTEXTTAU_mean = mean(TOTEXTTAU, na.rm = TRUE),
    T050_mean = mean(T050, na.rm = TRUE),
    TREFHT_mean = mean(TREFHT, na.rm = TRUE),
    SWTNT_mean = mean(SWTNT, na.rm = TRUE),
    LWTUP_mean = mean(LWTUP, na.rm = TRUE),
    SWGDNCLR_mean = mean(SWGDNCLR,na.rm=TRUE),
    TOTEXTTAU_sd = sd(TOTEXTTAU, na.rm = TRUE),
    T050_sd = sd(T050, na.rm = TRUE),
    TREFHT_sd = sd(TREFHT, na.rm = TRUE),
    SWTNT_sd = sd(SWTNT, na.rm = TRUE),
    LWTUP_sd = sd(LWTUP, na.rm = TRUE),
    SWGDNCLR_sd = sd(SWGDNCLR,na.rm=TRUE),
    .groups = "drop"
  )


# compute anomaleis
merra2_esn <-
  left_join(merra2, merra2_clmtg_spatial_means_and_sds, by = c("month", "lon", "lat")) %>%
  filter(abs(lat) <= 67.5) %>% # remove AEROD_v NAs
  mutate(
    TOTEXTTAU_stndzd = (TOTEXTTAU - TOTEXTTAU_mean) / (TOTEXTTAU_sd),
    LWTUP_stndzd = (LWTUP - LWTUP_mean) / (LWTUP_sd),
    SWTNT_stndzd = (SWTNT - SWTNT_mean),# / (SWTNT_sd),
    SWGDNCLR_stndzd = (SWGDNCLR - SWGDNCLR_mean) / (SWGDNCLR_sd),
    T050_stndzd = (T050 - T050_mean) / (T050_sd),
    TREFHT_stndzd = (TREFHT - TREFHT_mean) / (TREFHT_sd)
  )






# get order of lat/lon for weights
a <- merra2_esn %>%
  filter(data=="train") %>%
  select(date, lon, lat, TREFHT_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = TREFHT_stndzd
  ) 

lat_flat <- a$lat
weights <- cos(pi/180*lat_flat)

# create training/test matrices
train_mat_TOTEXTTAU <-
  merra2_esn %>%
  filter(data=="train") %>%
  select(date, lon, lat, TOTEXTTAU_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = TOTEXTTAU_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()



train_mat_LWTUP <-
  merra2_esn %>%
  filter(data=="train") %>%
  select(date, lon, lat, LWTUP_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = LWTUP_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_SWTNT <-
  merra2_esn %>%
  filter(data=="train") %>%
  select(date, lon, lat, SWTNT_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = SWTNT_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_SWGDNCLR <-
  merra2_esn %>%
  select(date, lon, lat, SWGDNCLR_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = SWGDNCLR_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_T050 <-
  merra2_esn %>%
  filter(data=="train") %>%
  select(date, lon, lat, T050_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = T050_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_TREFHT <-
  merra2_esn %>%
  filter(data=="train") %>%
  select(date, lon, lat, TREFHT_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = TREFHT_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

#
test_mat_TOTEXTTAU <-
  merra2_esn %>%
  filter(data=="test") %>%
  select(date, lon, lat, TOTEXTTAU_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = TOTEXTTAU_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()



test_mat_LWTUP <-
  merra2_esn %>%
  filter(data=="test") %>%
  select(date, lon, lat, LWTUP_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = LWTUP_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

test_mat_SWTNT <-
  merra2_esn %>%
  filter(data=="test") %>%
  select(date, lon, lat, SWTNT_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = SWTNT_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

test_mat_SWGDNCLR <-
  merra2_esn %>%
  filter(data=="test") %>%
  select(date, lon, lat, SWGDNCLR_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = SWGDNCLR_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

test_mat_T050 <-
  merra2_esn %>%
  filter(data=="test") %>%
  select(date, lon, lat, T050_stndzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = T050_stndzd
  ) %>%
  select(-lon,-lat) %>%
  t()

test_mat_TREFHT <-
  merra2_esn %>%
  filter(data=="test") %>%
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
eofs_TOTEXTTAU = compute_eofs(Ztrain = train_mat_TOTEXTTAU, n_eofs = n_eofs)
eofs_SWTNT = compute_eofs(Ztrain = train_mat_SWTNT, n_eofs = n_eofs)
eofs_LWTUP = compute_eofs(Ztrain = train_mat_LWTUP, n_eofs = n_eofs)
eofs_T050 = compute_eofs(Ztrain = train_mat_T050, n_eofs = n_eofs)
eofs_TREFHT = compute_eofs(Ztrain = train_mat_TREFHT, n_eofs = n_eofs)



# Model matrices
x_T050 = cbind(eofs_TOTEXTTAU$train, eofs_LWTUP$train, eofs_T050$train)
y_T050 = eofs_T050$train
x_T050_te = cbind(test_mat_TOTEXTTAU %*% eofs_TOTEXTTAU$phi, test_mat_LWTUP %*% eofs_LWTUP$phi, test_mat_T050 %*% eofs_T050$phi)

x_TREFHT = cbind(eofs_TOTEXTTAU$train, eofs_SWTNT$train, eofs_TREFHT$train)
y_TREFHT = eofs_TREFHT$train
x_TREFHT_te = cbind(test_mat_TOTEXTTAU %*% eofs_TOTEXTTAU$phi, test_mat_SWTNT %*% eofs_SWTNT$phi, test_mat_TREFHT %*% eofs_TREFHT$phi)

# for convenience
y_train_mat_T050 <- train_mat_T050 
phi_train_T050 <- eofs_T050$phi 

y_train_mat_TREFHT <- train_mat_TREFHT
phi_train_TREFHT <- eofs_TREFHT$phi 

t_T050 = rownames(x_T050)
t_TREFHT = rownames(x_TREFHT)
t_T050_te = rownames(x_T050_te)
t_TREFHT_te = rownames(x_TREFHT_te)

# hyperparameter values
tuning_params <-
  list(
    nensm = 5,
    cores = 20,
    tau = 1,
    m = 3,
    tau_emb = 1,
    nh = c(50, 100,200),
    U_width = c(0.1, 0.5),
    W_width = c(0.1, 0.5),
    U_pi = c(0.1, 0.5),
    W_pi = c(0.1, 0.5),
    nu = c(0.1, 0.5),
    reg_par = c(5,50,500, 1000),
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

write.csv(hs_pcs10_T050$output, file = "../output/gmd_paper/101-MERRA2_T050_hs_20pcs.csv", row.names = FALSE)


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

write.csv(hs_pcs10_TREFHT$output, file = "../output/gmd_paper/101-MERRA2_TREFHT_hs_20pcs.csv", row.names = FALSE)


