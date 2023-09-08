library(tidyverse)
library(rstan)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)


## Data loading

targets <- read_csv(here::here("data", "targets.csv"))


targets %>% 
  filter(Year == 2019) %>% data.frame()


drug <- targets %>% filter(Index == "DrugTime") %>% 
  filter(Year == 2019)

tx <- targets %>% filter(Index == "PrTxiPub") %>% 
  mutate(X = round(M * N))

det <- targets %>% 
  filter(Year == 2019) %>% 
  filter(Index == "TxI") %>% 
  mutate(
    N_Txi = round(N * M),
    N_Det = N
  )


ds <- local({
  N_Test_SSM = 13914911
  N_Test_Xpert = 4120552 
  N_Det_Pub = 1688427
  N_Det_Eng = 733694
  N_DetBac = 513050 + 548981
  N_DetCDx = 835930 + 1037210
  N_Txi_Pub = 1527464
  N_Txi_Eng = 554980
  N_Pop = 1425775850
  
  PropXpert_Eng = (215594 + 262160) / (3483130 + 2365739)
  # BacPerXpert_Eng = 68556 / (215594 + 262160)
  
  N_Test_Xpert_Eng = round(N_Test_Xpert * PropXpert_Eng)
  # N_DetBac_Eng = round(N_Test_Xpert_Eng * BacPerXpert_Eng)
  
  list(
    N_Det_Pub = N_Det_Pub,
    N_Det_Eng = N_Det_Eng,
    N_DetBac = N_DetBac,
    # N_DetBac_Pub = N_DetBac - N_DetBac_Eng,
    # N_DetBac_Eng = N_DetBac_Eng,
    N_Test_SSM_Pub = N_Test_SSM,
    N_Test_NAAT_Pub = N_Test_Xpert - N_Test_Xpert_Eng,
    N_Test_NAAT_Eng = N_Test_Xpert_Eng,
    N_DetCDx = N_Det_Pub + N_Det_Eng - N_DetBac,
    # N_DetCDx_Pub = N_Det_Pub - (N_DetBac - N_DetBac_Eng),
    # N_DetCDx_Eng = N_Det_Eng - N_DetBac_Eng,
    N_Txi_Pub = N_Txi_Pub,
    N_Txi_Eng = N_Txi_Eng,
    
    Pop = N_Pop,
    
    tp_bac = c(0.73984, 0.544, 0.7225, 0),
    fp_bac = c(0.0336, 0.017, 0.017, 0),
    test_tb_ssm = c(0.85, 0.85, 0, 0),
    test_tb_naat = c(0.306, 0, 0.85, 0),
    test_nontb_ssm = c(0.85, 0.85, 0, 0),
    test_nontb_naat = c(0.833, 0, 0.85, 0),
    
    Tx = tx$N,
    Tx_Pub = tx$X,
    Drug = drug$M,
    Drug_Std = drug$Error,
    
    p_csi_pub = 0.483,
    dur_upper = 2
  )
})

