library(tidyverse)
library(tidybayes)


theme_set(theme_bw())


scs <- c(
  Baseline = "Baseline",
  RedRel = "Reduce post-treatment relapse",
  Txi = "Reduce pre-treatment LTFU",
  TSwab = "Tongue swab with current Xpert",
  MassScreen = "Mass screening",
  Vac = "Post-exposure vaccination",
  All = "All together"
)



sim <- read_csv(here::here("out", "post_free", "Sim_Intv.csv")) %>% 
  filter(Year <= 2035 & Year >= 2024)




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
  )


g_epi <- sim %>% 
  select(Year, IncR, MorR, Scenario) %>% 
  pivot_longer(c(IncR, MorR)) %>% 
  group_by(Year, Scenario, name) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  mutate(
    Scenario = scs[Scenario],
    Scenario = factor(Scenario, scs)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
  geom_line(aes(x = Year, y = M, colour = Scenario)) +
  scale_y_continuous("per 100k", labels = scales::number_format(scale = 1e5)) +
  scale_x_continuous("Year", breaks = c(2024, seq(2025, 2035, 5))) + 
  facet_wrap(.~name, scale = "free_y", 
             labeller = labeller(name = c("IncR" = "Incidence", "MorR" = "Mortality"))) +
  expand_limits(y = 0)

g_epi


g_avt <- avt %>% 
  select(Year, AvtInc, AvtMor, Scenario) %>% 
  pivot_longer(c(AvtInc, AvtMor)) %>% 
  group_by(Year, Scenario, name) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  mutate(
    Scenario = scs[Scenario],
    Scenario = factor(Scenario, scs)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
  geom_line(aes(x = Year, y = M, colour = Scenario)) +
  scale_x_continuous("Year", breaks = c(2024, seq(2025, 2035, 5))) + 
  scale_y_continuous("per cent", labels = scales::percent) +
  facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
  expand_limits(y = 0)

g_avt


sim %>% 
  select(Year, IncRecentR, IncRemoteR, Scenario) %>% 
  pivot_longer(c(IncRecentR, IncRemoteR)) %>% 
  group_by(Year, Scenario, name) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  mutate(
    Scenario = scs[Scenario],
    Scenario = factor(Scenario, scs)
  ) %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
  geom_line(aes(x = Year, y = M, colour = Scenario)) +
  facet_wrap(.~name, scale = "free_y") +
  scale_x_continuous("Year", breaks = c(2024, seq(2025, 2035, 5))) + 
  scale_y_continuous("per 100k", labels = scales::number_format(scale = 1e5)) +
  expand_limits(y = 0)


g_avt

ggsave(g_epi, filename = here::here("docs", "g_epi.png"), width = 10, height = 5)
ggsave(g_avt, filename = here::here("docs", "g_avt.png"), width = 10, height = 5)


