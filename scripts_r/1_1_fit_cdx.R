library(jsonlite)
library(rstan)

model <- rstan::stan_model(here::here("stan", "cdx.stan"))

folder <- "pars_cdx"
dir.create(here::here("out", folder), showWarnings = F)


dat_itr <- jsonlite::read_json(here::here("data", "itr.json"))


dat <- with(dat_itr, {
  pr_pri <- (Test_CBNAAT$Detected_Private + Test_Truenat$Detected_Private) / (Test_CBNAAT$Detected + Test_Truenat$Detected)
  
  d <- list(
    Pop = 1375586000,
    N_Pre_SSM = Presumptive$SSM,
    N_Pre_Xpert = Presumptive$CBNAAT + Presumptive$Truenat, 
    N_Det_SSM_Pub = Notification$Test_SSM,
    N_Det_Xpert_Pub = round(Notification$Test_Xpert * (1 - pr_pri)),
    N_Det_Xpert_Eng = round(Notification$Test_Xpert * pr_pri)
  )
  d$N_Det_CDX_Pub <- Notification$Sector_Public - d$N_Det_SSM_Pub - d$N_Det_Xpert_Pub
  d$N_Det_CDX_Eng <- Notification$Sector_Private - d$N_Det_Xpert_Eng
  
  d
})

save(dat, file = here::here("out", folder, "data_src.rdata"))



## Model fitting -----
post <- rstan::sampling(model, data=dat, iter=1e4, warmup=1e4 - 1000)

post

pars <- extract(post) %>% data.frame() %>% as.tibble()

write_csv(pars, file = here::here("pars", "pars_cdx.csv"))

write_json(pars, path = here::here("pars", "pars_cdx.json"))


