library(tidyverse)
library(rstan)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)


## Exogenous variables

exo <- jsonlite::read_json(here::here("data", "pre_cdx.json"), simplifyVector=T)


## Data loading

targets <- read_csv(here::here("data", "targets.csv"))


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
  
  
   det <- targets %>% 
     filter(Index == "CNR") %>% 
     filter(Year %in% c(2021, 2022)) %>% 
     filter(Tag != "All") %>% 
     mutate(
       Tag = case_when(
         Tag == "NAAT" ~ "Xpert",
         Tag == "Pri" ~ "Eng",
         T~ Tag
       ),
       Key = paste0("N_Det_", Tag),
       X = round(N * M)
     ) %>% 
     select(Key, Year, X) %>% 
     pivot_wider(values_from = X, names_from = Key) %>% 
     arrange(Year) %>% 
     select(- Year) %>% 
     as.list()
  
  
  res <- c(res, det)
  
  
  tests <- targets %>% 
    filter(Index == "TestR") %>% 
    filter(Year %in% c(2021, 2022)) %>% 
    filter(Tag != "All") %>% 
    mutate(
      Key = paste0("N_Test_", Tag),
      X = round(N * M)
    ) %>% 
    select(Key, Year, X) %>% 
    arrange(Year) %>% 
    pivot_wider(values_from = X, names_from = Key) %>% 
    arrange(Year) %>% 
    mutate(
      PrXpertPri = c((155295 + 228783) / (2197757 + 1434124), 
                     (215594 + 262160) / (3483130 + 2365739)),
      N_Test_Xpert_Pub = round(N_Test_NAAT * (1 - PrXpertPri)),
      N_Test_Xpert_Eng = round(N_Test_NAAT * PrXpertPri)
    ) %>% 
    select(N_Test_SSM_Pub = N_Test_SSM, N_Test_Xpert_Pub, N_Test_Xpert_Eng) %>% 
    as.list() 

  res <- c(res, tests)
  
  ## Todo
  txi <- list(
    N_Txi_Pub = c(1275823, 1527464),
    N_Txi_Eng = c(475614, 554980)
  )
  res <- c(res, txi)
  
  res <- c(res, exo)
  res
})


## Extract target data2plot

targets2plot <- bind_rows(
    as_tibble(ds[startsWith(names(ds), "N_")]) %>% 
      mutate(
        Year = ds$Years,
        Pop = ds$Pop
      ) %>% 
      pivot_longer(-c(Year, Pop)) %>% 
      mutate(
        m = value / Pop,
        key = gsub("N_", "", name)
      ) %>% 
      select(Year, key, m, value),
    tribble(
      ~Year, ~key, ~m, ~value,
      2021, "Drugtime", ds$Drug, ds$Drug,
    )
  )






tar <- read_csv(here::here("data", "targets_india.csv"))

prev <- local({
  prev <- targets %>% filter(startsWith(Index, "Prev") | Index == "PrCSIPub" | Index == "TBLikeUt")
  prev <- as.list(setNames(prev$M, prev$Index))
})


txo <- targets %>% 
  filter(Index == "TxSucc" | Index == "TxDie") %>% 
  group_by(Index, Tag) %>% 
  summarise(M = weighted.mean(M, N))



## Model fitting

model_key = "cs_independent"

for(model_key in c("cs_free")) {
  
  model <- rstan::stan_model(here::here("scripts", glue::as_glue(model_key) + ".stan"))
  
  post <- rstan::sampling(model, data=ds, iter=1e4, warmup=1e4 - 1000)
  
  stan_dens(post, pars=c("p_tb_pub", "p_tb_eng", "p_tb_pri", "r_test"))
  
  stan_dens(post, pars=c("p_ava_ssm_pub", "p_ava_xpert_pub", "p_ava_xpert_eng"))
  stan_dens(post, pars=c("alg_pub"))
  pairs(post, pars=c("alg_pub"))
  stan_dens(post, pars=c("alg_eng"))
  
  stan_dens(post, pars=c("ppv_pub", "ppv_eng", "ppv_pri"))
  stan_scat(post, pars=c("sens_cdx", "spec_cdx"))
  
  
  stan_scat(post, pars=c("dur_pri", "p_txi_pri"))
  
  stan_dens(post, pars=c("p_csi_ppm"))
  
  
  ext <- extract(post) %>% 
    data.frame() %>% 
    as.tibble() %>% 
    mutate(
      p_csi_pub = ds$p_csi_pub
    )
  
  tab <- as.data.frame(summary(post)$summary)
  tab$Name <- rownames(tab)
  tab <- tab %>% as_tibble() %>% relocate(Name)
  
  
  dir.create(here::here("out", model_key), showWarnings = F)
  save(post, file = here::here("out", model_key, "post.rdata"))
  write_csv(tab, file = here::here("out", model_key, "summary.csv"))
  write_csv(ext, file = here::here("out", model_key, "post.csv"))
  
  
  to_fit <- ext %>% 
    mutate(
      R_DetPub = r_det_bac_pub + r_det_cdx_pub,
      R_DetEng = r_det_bac_eng + r_det_cdx_eng,
      R_DrugTime = drug_time,
      R_TestSSM = r_test_ssm_pub,
      R_TestNAAT = r_test_naat_pub + r_test_naat_eng,
      R_NotiBac = r_det_bac_pub + r_det_bac_eng,
      R_NotiCDx = r_det_cdx_pub + r_det_cdx_eng,
      R_TxiPub = (r_det_bac_pub + r_det_cdx_pub) * p_txi_pub,
      R_TxiEng = (r_det_bac_eng + r_det_cdx_eng) * p_txi_eng
    ) %>% 
    select(starts_with("R", ignore.case = F)) %>% 
    pivot_longer(everything(), names_to = "Index") %>% 
    mutate(
      Index = gsub("R_", "", Index)
    )
  
  
  
  g_gof <- to_fit %>% 
    ggplot() +
    geom_density(aes(x = value)) +
    geom_vline(data = targets, aes(xintercept = value)) + 
    facet_wrap(.~Index, scales = "free_x") +
    scale_x_continuous("rate, per 100k", labels = scales::number_format(scale = 1e5))
  
  g_gof
  
  ggsave(g_gof, filename = here::here("out", model_key, "g_gof.png"), width = 7, height = 6)
  
  
  sel <- ext[sample.int(nrow(ext), 2000), ]
  
  js <- list(
    pars = sel,
    prev = prev
  )
  
  jsonlite::write_json(js, here::here("docs", "pars", "pars_" + glue::as_glue(model_key) + ".json"), digits = 8, auto_unbox = T)
  
}



ext <- extract(post) %>% 
  data.frame() %>% 
  as.tibble() %>% 
  mutate(
    p_csi_pub = ds$p_csi_pub
  ) %>% 
  select(
    starts_with(c("r_det_pub", "r_det_eng", "r_det_bac.", "r_det_cdx.", "r_test_ssm", "r_test_xpert", "r_txi_pub", "r_txi_eng")),
    r_DrugTime.1 = drug_time
  )


to_fit <- ext %>% 
  mutate(
    k = 1:n()
  ) %>% 
  filter(k <= 500) %>% 
  pivot_longer(-k, 
               names_to = c("Index", "Year"),
               names_pattern = "r_(\\w+).(1|2)") %>% 
  mutate(
    Index = stringr::str_to_title(Index),
    Index = gsub("ssm", "SSM", Index),
    Index = gsub("pub", "Pub", Index),
    Index = gsub("eng", "Eng", Index),
    Index = gsub("xpert", "Xpert", Index),
    Index = gsub("bac", "Bac", Index),
    Index = gsub("cdx", "CDx", Index),
    Year = ifelse(Year == 1, 2021, 2022)
  )


to_fit %>% 
  group_by(Index, Year) %>% 
  summarise(
    m = median(value),
    l = quantile(value, 0.1),
    u = quantile(value, 0.9)
  ) %>% 
  left_join(targets2plot %>% select(Index = key, Year, d = m))


to_fit %>% 
  ggplot() +
  geom_density(aes(x = value, fill = as.factor(Year))) +
  geom_vline(data = targets2plot %>% select(Index = key, Year, m), aes(xintercept = m)) + 
  facet_wrap(Index~Year, scales = "free_x", ncol = 4)



