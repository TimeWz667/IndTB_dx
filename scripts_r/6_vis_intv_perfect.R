library(tidyverse)
library(tidybayes)


theme_set(theme_bw())


qci <- 0.5


## Loading data -----
scs <- c(
  Baseline = "Baseline",
  PerfectDx = "No Misdiagnosis",
  PerfectTxi = "No Pretreatment LTFU",
  NoRelapse = "No Relapse after treatment"
)

sim <- read_csv(here::here("out", "post_dyage", "Sim_IntvPerfect.csv"))
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
  )


g_avt <- avt %>% 
  select(Year, AvtInc, AvtMor, Scenario) %>% 
  pivot_longer(c(AvtInc, AvtMor)) %>% 
  group_by(Year, Scenario, name) %>% 
  summarise(
    M = median(value),
    L = quantile(value, (1 - qci) / 2),
    U = quantile(value, 1 - (1 - qci) / 2)
  ) %>% 
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

g_avt
ggsave(g_avt, filename = here::here("docs", "figs", "g_intv_perfect_avt.png"), width = 12, height = 5)

