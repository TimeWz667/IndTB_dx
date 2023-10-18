library(jsonlite)
library(rstan)

model <- rstan::stan_model(here::here("stan", "dx.stan"))

folder <- "pars_dx"
dir.create(here::here("out", folder), showWarnings = F)


pars_dx <- jsonlite::read_json(here::here("data", "dx.json"))
dat_itr <- jsonlite::read_json(here::here("data", "itr.json"))


dat <- c(
  with(dat_itr, {
    pr_pri <- (Test_CBNAAT$Tested_Private + Test_Truenat$Tested_Private) / (Test_CBNAAT$Tested + Test_Truenat$Tested)
    
    list(
      N_Pre_SSM = Presumptive$SSM,
      N_Pre_Xpert = Presumptive$CBNAAT + Presumptive$Truenat, 
      N_Det_SSM = Notification$Test_SSM,
      N_Det_Xpert = Notification$Test_Xpert,
      N_Det_CDX = Notification$Test_CXR + Notification$Test_Other,
      N_Det_Pub = Notification$Sector_Public,
      N_Det_Eng = Notification$Sector_Private,
      PrEng = pr_pri
    )
  }),
  with(pars_dx, {
    list(
      Sens_SSM = sens_ssm,
      Sens_Xpert = sens_xpert,
      Sens_XpertSN = pars_dx[["sens_xpert_ss-"]],
      Spec_SSM = spec_ssm,
      Spec_Xpert = spec_xpert,
      Spec_XpertSN = spec_xpert
    )
  })
)


dat$Pop <- 1375586000
save(dat, file = here::here("out", folder, "data_src.rdata"))



## Model fitting -----
post <- rstan::sampling(model, data=dat, iter=1e4, warmup=1e4 - 1000)

post


