library(tidyverse)
library(rstan)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)


## Model loading

model <- rstan::stan_model(here::here("stan", "bac_cdx_cs.stan"))


## Data loading

targets <- read_csv(here::here("data", "targets.csv"))

prev <- local({
  prev <- targets %>% filter(startsWith(Index, "Prev") | Index == "PrCSIPub" | Index == "TBLikeUt")
  prev <- as.list(setNames(prev$M, prev$Index))
})


txo <- targets %>% 
  filter(Index == "TxSucc" | Index == "TxDie") %>% 
  group_by(Index, Tag) %>% 
  summarise(M = weighted.mean(M, N))


txi <- list(
  "bac_cdx_cs_2021" = list(
    N_Txi_Pub = 1275823, N_Txi_Eng = 475614
  ),
  "bac_cdx_cs_2022" = list(
    N_Txi_Pub = 1527464, N_Txi_Eng = 554980
  )
)


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
  )

## 

for (yr in c(2021, 2022)) {
  folder <- "bac_cdx_cs_" + glue::as_glue(yr)
  dir.create(here::here("out", folder), showWarnings = F)
  
  ds <- local({
    tx <- targets %>% filter(Index == "PrTxiPub") %>% 
      mutate(X = round(M * N))
    
    drug <- targets %>% filter(Index == "DrugTime") %>% 
      filter(Year == 2019)
    
    n_pop <- targets %>% 
      filter(Index == "CNR" & Tag == "All") %>% 
      filter(Year == yr) %>% pull(N)
    
    res <- list(
      Pop = n_pop,
      Tx = tx$N,
      Tx_Pub = tx$X,
      Drug = drug$M,
      Drug_Std = drug$Error,
      # p_csi_pub = 0.483,
      dur_upper = 2
    )
    
    
    itr <- c(
      targets %>% 
        filter(Index == "CNR") %>% 
        filter(Year == yr) %>% 
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
        select(Key, X) %>% 
        pivot_wider(values_from = X, names_from = Key) %>% 
        as.list(),
        tests %>%
        filter(Year == yr) %>% 
        select(N_Test_SSM_Pub = N_Test_SSM, N_Test_Xpert_Pub, N_Test_Xpert_Eng) %>% 
        as.list()
    )
    
    ds <- c(res, itr, exo)
    ds <- c(ds, txi[[folder]])
    ds
  })
  
  
  save(ds, file = here::here("out", folder, "data_src.rdata"))
  
  
  ## Model fitting
  post <- rstan::sampling(model, data=ds, iter=1e4, warmup=1e4 - 1000)
  
  ext <- extract(post) %>% 
    data.frame() %>% 
    as.tibble()
  
  
  ## GOF
  stan_dens(post, pars=c("alg_pub", "alg_eng"))
  stan_scat(post, pars=c("sens_cdx", "spec_cdx"))
  
  targets2plot <- bind_rows(
    as_tibble(ds[startsWith(names(ds), "N_")]) %>% 
      mutate(
        Pop = ds$Pop
      ) %>% 
      pivot_longer(-c(Pop)) %>% 
      mutate(
        m = value / Pop,
        Index = gsub("N_", "", name)
      ) %>% 
      select(Index, m, value),
    tribble(
      ~Index, ~m, ~value,
      "Drugtime", ds$Drug, ds$Drug,
    )
  )
  
  
  to_fit <- ext %>% 
    mutate(
      R_Det_Pub = r_det_pub,
      R_Det_Eng = r_det_eng,
      R_Test_SSM_Pub = r_test_ssm_pub,
      R_Test_Xpert_Pub = r_test_xpert_pub,
      R_Test_Xpert_Eng = r_test_xpert_eng,
      R_Det_Bac = r_det_bac,
      R_Det_CDx = r_det_cdx,
      R_Txi_Pub = (r_det_pub) * p_txi_pub,
      R_Txi_Eng = (r_det_eng) * p_txi_eng,
      R_Drugtime = drug_time
    ) %>% 
    select(starts_with("R", ignore.case = F)) %>% 
    pivot_longer(everything(), names_to = "Index") %>% 
    mutate(
      Index = gsub("R_", "", Index)
    )
  
  g_gof <- to_fit %>% 
    ggplot() +
    geom_density(aes(x = value)) +
    geom_vline(data = targets2plot, aes(xintercept = m)) + 
    facet_wrap(.~Index, scales = "free_x") +
    scale_x_continuous("rate, per 100k", labels = scales::number_format(scale = 1e5)) +
    labs(subtitle = yr)
  
  g_gof
  ggsave(g_gof, filename = here::here("out", folder, "g_gof.png"), width = 7, height = 6)
  
  
  ## Output pars
  tab <- as.data.frame(summary(post)$summary)
  tab$Name <- rownames(tab)
  tab <- tab %>% as_tibble() %>% relocate(Name)
  
  
  save(post, file = here::here("out", folder, "post.rdata"))
  write_csv(tab, file = here::here("out", folder, "summary.csv"))
  write_csv(ext, file = here::here("out", folder, "post.csv"))
  
  
  ## To Parameters
  sel <- ext[sample.int(nrow(ext), 2000), ]
  
  js <- list(
    pars = sel,
    prev = prev,
    txo = txo
  )
  
  jsonlite::write_json(js, here::here("pars", "pars_" + folder + ".json"), digits = 8, auto_unbox = T)

}


