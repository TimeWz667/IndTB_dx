library(tidyverse)
library(rstan)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)


## Exogenous variables

exo <- jsonlite::read_json(here::here("data", "pre_cdx.json"), simplifyVector=T)


## Data loading

targets <- read_csv(here::here("data", "targets.csv"))


det <- targets %>% 
  filter(Year == 2019) %>% 
  filter(Index == "TxI") %>% 
  mutate(
    N_Txi = round(N * M),
    N_Det = N
  )


ds <- local({
  n_pop <- targets %>% 
    filter(Index == "CNR" & Tag == "All") %>% 
    filter(Year %in% c(2021, 2022)) %>% pull(N)
  
  tx <- targets %>% filter(Index == "PrTxiPub") %>% 
    mutate(X = round(M * N))
  
  
  drug <- targets %>% filter(Index == "DrugTime") %>% 
    filter(Year == 2019)
  
  res <- list(
    Year0 = 2020,
    Years = c(2021, 2022),
    Pop = n_pop,
    Tx = tx$N,
    Tx_Pub = tx$X,
    Drug = drug$M,
    Drug_Std = drug$Error,
    p_csi_pub = 0.483,
    dur_upper = 2
  )
  
  
  det <- local({
    temp <- targets %>% 
      filter(Index == "CNR") %>% 
      filter(Year %in% c(2021, 2022)) %>% 
      filter(Tag != "All") %>% 
      mutate(
        Key = paste0("N_Det_", Tag),
        X = round(N * M)
      ) %>% 
      select(Key, Year, X) %>% 
      arrange(Year)
    
    temp
    
    ks <- unique(temp$Key)
    ks
    
    lapply(set_names(ks, ks), function(k) temp %>% filter(Key == k) %>% pull(X))
  })
  
  
  res <- c(res, det)
  
  
  tests <- local({
    temp <- targets %>% 
      filter(Index == "TestR") %>% 
      filter(Year %in% c(2021, 2022)) %>% 
      filter(Tag != "All") %>% 
      mutate(
        Key = paste0("N_Test_", Tag),
        X = round(N * M)
      ) %>% 
      select(Key, Year, X) %>% 
      arrange(Year)
    
    temp
    
    ks <- unique(temp$Key)
    ks
    
    lapply(set_names(ks, ks), function(k) temp %>% filter(Key == k) %>% pull(X))
  })
  
  
  ds <- list(
    N_Det_Pub = N_Det_Pub,
    N_Det_Eng = N_Det_Eng,
    N_DetBac = N_DetBac,
    # N_DetBac_Pub = N_DetBac - N_DetBac_Eng,
    # N_DetBac_Eng = N_DetBac_Eng,
    N_Test_SSM_Pub = N_Test_SSM,
    N_Test_Xpert_Pub = N_Test_Xpert - N_Test_Xpert_Eng,
    N_Test_Xpert_Eng = N_Test_Xpert_Eng,
    N_DetCDx = N_Det_Pub + N_Det_Eng - N_DetBac,
    # N_DetCDx_Pub = N_Det_Pub - (N_DetBac - N_DetBac_Eng),
    # N_DetCDx_Eng = N_Det_Eng - N_DetBac_Eng,
    N_Txi_Pub = N_Txi_Pub,
    N_Txi_Eng = N_Txi_Eng,
    

  )
  
  res <- c(res, exo)
  res
})

