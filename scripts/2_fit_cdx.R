library(tidyverse)
library(rstan)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)


## Model loading

model <- rstan::stan_model(here::here("stan", "bac_cdx.stan"))



## Exogenous variables

exo <- jsonlite::read_json(here::here("data", "pre_cdx.json"), simplifyVector=T)


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



for (yr in c(2021, 2022)) {
  folder <- "bac_cdx_" + glue::as_glue(yr)
  dir.create(here::here("out", folder), showWarnings = F)

  ds <- c(
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
      select(Key, Year, X) %>% 
      pivot_wider(values_from = X, names_from = Key) %>% 
      select(N_Det_Bac, N_Det_CDx, N_Det_Xpert) %>% 
      as.list(),
    targets %>% 
      filter(Index == "TestR") %>% 
      filter(Year == yr) %>% 
      filter(Tag != "All") %>% 
      mutate(
        Key = paste0("N_Test_", Tag),
        X = round(N * M)
      ) %>% 
      select(Key, X) %>% 
      pivot_wider(values_from = X, names_from = Key) %>% 
      select(N_Test_SSM, N_Test_Xpert = N_Test_NAAT) %>% 
      as.list()
  )
  
  
  ds$Pop <- targets %>% 
    filter(Index == "CNR" & Tag == "All") %>% 
    filter(Year == yr) %>% pull(N)
  
  ds <- c(ds, exo)
  save(ds, file = here::here("out", folder, "data_src.rdata"))
  
  ## Model fitting
  
  post <- rstan::sampling(model, data=ds, iter=1e4, warmup=1e4 - 1000)
  
  ext <- extract(post) %>% 
    data.frame() %>% 
    as.tibble() %>% 
    mutate(
      p_csi_pub = ds$p_csi_pub
    )
  
  ## GOF
  
  stan_dens(post, pars=c("alg"))
  stan_scat(post, pars=c("sens_cdx", "spec_cdx"))
  
  targets2plot <- as_tibble(ds[startsWith(names(ds), "N_")]) %>% 
    mutate(
      Pop = ds$Pop
    ) %>% 
    pivot_longer(-c(Pop)) %>% 
    mutate(
      m = value / Pop,
      Index = gsub("N_", "", name)
    ) %>% 
    select(Index, m, value)
  
  
  to_fit <- ext %>% 
    mutate(
      R_Test_SSM = r_test_ssm,
      R_Test_Xpert = r_test_xpert,
      R_Det_Bac = r_det_bac,
      R_Det_CDx = r_det_cdx
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
  
  ## Posterior output
  
  
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



