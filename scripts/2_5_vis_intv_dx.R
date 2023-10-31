library(tidyverse)
library(tidybayes)


theme_set(theme_bw())


qci <- 0.5


scs <- c(
  Baseline = "Baseline",
  Dx_TSwab = "Swab with current tests",
  Dx_POC = "Swab with PoC test",
  Dx_POC_Hi = "Swab with\nhigh-throughput NAAT"
)

## Loading data -----

for (folder in c("post_dy_c0", "post_dy_ca")) {
  sim <- read_csv(here::here("out", folder, "Sim_Intv.csv"))
  sim <- sim %>% filter(Year >= 2024 & Year <= 2040)
  
  
  tab_epi <- sim %>% 
    select(Year, IncR, MorR, Scenario) %>% 
    pivot_longer(c(IncR, MorR)) %>% 
    group_by(Year, Scenario, name) %>% 
    summarise(
      M = median(value),
      L = quantile(value, (1 - qci) / 2),
      U = quantile(value, 1 - (1 - qci) / 2)
    )
  
  g_epi <- tab_epi %>%
    filter(Scenario %in% names(scs)) %>% 
    mutate(
      PPM = ifelse(str_detect(Scenario, "PPM"), "With PPM", "Without PPM"),
      Scenario = gsub("PPM", "", Scenario)
    ) %>% 
    mutate(
      Scenario = scs[Scenario],
      Scenario = factor(Scenario, scs)
    ) %>% 
    ggplot() +
    geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
    geom_line(aes(x = Year, y = M, colour = Scenario)) +
    scale_y_continuous("per 100k", labels = scales::number_format(scale = 1e5)) +
    scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
    facet_wrap(.~name, scale = "free_y", 
               labeller = labeller(name = c("IncR" = "Incidence", "MorR" = "Mortality"))) +
    expand_limits(y = 0)
  
  
  tab_avt <- sim %>% 
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
  
  
  g_avt <- tab_avt %>% 
    filter(Scenario %in% names(scs)) %>% 
    mutate(
      Scenario = scs[Scenario],
      Scenario = factor(Scenario, scs)
    ) %>% 
    ggplot() +
    geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
    geom_line(aes(x = Year, y = M, colour = Scenario)) +
    scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
    scale_y_continuous("per cent", labels = scales::percent) +
    facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
    expand_limits(y = 0)
  
  ggsave(g_epi, filename = here::here("out", folder, "g_intv_dx_epi.png"), width = 12, height = 5)
  ggsave(g_avt, filename = here::here("out", folder, "g_intv_dx_avt.png"), width = 12, height = 5)
  
  
}



tab <- bind_rows(lapply(c("post_dy_c0", "post_dy_ca"), function(folder) {
  sim <- read_csv(here::here("out", folder, "Sim_Intv.csv"))
  sim <- sim %>% filter(Year >= 2024 & Year <= 2040)
  
  tab_avt <- sim %>% 
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
  
  tab_avt %>% mutate(Assumption = folder) %>% 
    filter(Scenario %in% names(scs)) %>% 
    mutate(
      Scenario = scs[Scenario],
      Scenario = factor(Scenario, scs)
    )
}))


g_avtinc <- tab %>% 
  mutate(
    Assumption = ifelse(Assumption == "post_dy_c0", "1CDx", "2CDx"),
    Assumption = factor(Assumption, c("2CDx", "1CDx"))
  ) %>% 
  filter(name == "AvtInc") %>% 
  ggplot() +
  geom_line(aes(x = Year, y = M, colour = Scenario, linetype=Assumption)) +
  scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
  scale_y_continuous("Averted incidence, %", labels = scales::percent) +
  #facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
  expand_limits(y = 0)

g_avtmor <- tab %>% 
  mutate(
    Assumption = ifelse(Assumption == "post_dy_c0", "1CDx", "2CDx"),
    Assumption = factor(Assumption, c("2CDx", "1CDx"))
  ) %>% 
  filter(name == "AvtMor") %>% 
  ggplot() +
  geom_line(aes(x = Year, y = M, colour = Scenario, linetype=Assumption)) +
  scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
  scale_y_continuous("Averted mortality, %", labels = scales::percent) +
  #facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
  expand_limits(y = 0)


ggsave(g_avtinc, filename = here::here("docs", "figs", "g_intv_dx_avtinc.png"), width = 6, height = 4)
ggsave(g_avtmor, filename = here::here("docs", "figs", "g_intv_dx_avtmor.png"), width = 6, height = 4)

