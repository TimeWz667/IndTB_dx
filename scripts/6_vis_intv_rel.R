library(tidyverse)
library(tidybayes)


theme_set(theme_bw())


qci <- 0.5


## Loading data -----

sim <- read_csv(here::here("out", "post_dyage", "Sim_IntvRel.csv"))
sim <- sim %>% filter(Year >= 2024 & Year <= 2040)


avt <- sim %>% 
  select(Year, CumInc, CumMor, Scenario, Key) %>% 
  left_join(
    sim  %>% 
      filter(Scenario == "Baseline") %>% 
      select(Year, CumInc0 = CumInc, CumMor0 = CumMor, Key)
  ) %>% 
  group_by(Key, Scenario) %>% 
  mutate(
    CumInc = CumInc - CumInc[1],
    CumMor = CumMor - CumMor[1],
    CumInc0 = CumInc0 - CumInc0[1],
    CumMor0 = CumMor0 - CumMor0[1]
  ) %>% 
  ungroup() %>% 
  mutate(
    DiffInc = CumInc0 - CumInc,
    DiffMor = CumMor0 - CumMor,
    AvtInc = DiffInc / CumInc0,
    AvtMor = DiffMor / CumMor0,
    AvtInc = ifelse(is.na(AvtInc), 0, AvtInc),
    AvtMor = ifelse(is.na(AvtMor), 0, AvtMor)
  ) %>% 
  select(Year, AvtInc, AvtMor, Scenario) %>% 
  pivot_longer(c(AvtInc, AvtMor)) %>% 
  group_by(Year, Scenario, name) %>% 
  summarise(
    M = median(value),
    L = quantile(value, (1 - qci) / 2),
    U = quantile(value, 1 - (1 - qci) / 2)
  ) 
