library(tidyverse)
library(tidybayes)


theme_set(theme_bw())


qci <- 0.5


scs <- c(
  Baseline = "Baseline",
  "Combine_Lo" = "Combination (conservative)",
  "Combine_Hi" = "Combination (ambitious)"
)



sim2 <- read_csv(here::here("out", "post_dyage", "Sim_IntvDelayed.csv")) %>% 
  filter(Year >= 2024 & Year <= 2040) %>% 
  filter(Scenario %in% c("Baseline", "Combine_Lo", "Combine_Hi")) %>% 
  mutate(Delayed = "2X")

sim <- read_csv(here::here("out", "post_dyage", "Sim_IntvAll.csv")) %>% 
  filter(Year >= 2024 & Year <= 2040) %>% 
  filter(Scenario %in% c("Baseline", "Combine_Lo", "Combine_Hi")) %>% 
  mutate(Delayed = "None")



avt <- sim %>% 
  select(Year, CumInc, CumMor, Scenario, Key, Delayed) %>% 
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
  ) %>% 
  filter(Scenario %in% names(scs)) %>% 
  mutate(
    Scenario = scs[Scenario],
    Scenario = factor(Scenario, scs)
  )

avt2 <- sim2 %>% 
  select(Year, CumInc, CumMor, Scenario, Key, Delayed) %>% 
  left_join(
    sim2  %>% 
      filter(Scenario == "Baseline") %>% 
      select(Year, CumInc0 = CumInc, CumMor0 = CumMor, Key)
  ) %>% 
  group_by(Key, Scenario, Delayed) %>% 
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
  select(Year, AvtInc, AvtMor, Scenario, Delayed) %>% 
  pivot_longer(c(AvtInc, AvtMor)) %>% 
  group_by(Year, Scenario, name, Delayed) %>% 
  summarise(
    M = median(value),
    L = quantile(value, (1 - qci) / 2),
    U = quantile(value, 1 - (1 - qci) / 2)
  ) %>% 
  filter(Scenario %in% names(scs)) %>% 
  mutate(
    Scenario = scs[Scenario],
    Scenario = factor(Scenario, scs)
  )


g_avt <- avt %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
  geom_line(aes(x = Year, y = M, colour = Scenario, linetype = "2 years")) +
  geom_line(data = avt2, aes(x = Year, y = M, colour = Scenario, linetype = "4 years"), linewidth = rel(1.2)) +
  geom_line(data = avt2, aes(x = Year, y = L, colour = Scenario, linetype = "4 years")) +
  geom_line(data = avt2, aes(x = Year, y = U, colour = Scenario, linetype = "4 years")) +
  scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
  scale_y_continuous("per cent", labels = scales::percent) +
  scale_linetype_discrete("Scale-up period") +
  facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
  expand_limits(y = 0)

g_avt
ggsave(g_avt, filename = here::here("docs", "figs", "g_intv_delayed_avt.png"), width = 12, height = 5)



