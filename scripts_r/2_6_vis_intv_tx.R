library(tidyverse)
library(tidybayes)


theme_set(theme_bw())


qci <- 0.5


scs <- c(
  Baseline = "Baseline",
  "Tx_BPaLM" = "BPaLM",
  "Tx_PAN-TB" = "PAN-TB",
  "Tx_LA-INJ" = "LA-INJ"
)

folder <- "dy"

sim <- read_csv(here::here("out", folder, "Sim_IntvTx.csv"))
sim <- sim %>% filter(Year >= 2024 & Year <= 2040)

sim %>% 
  group_by(Scenario, Year) %>% 
  summarise(V = mean(IncR)) %>% 
  ggplot() +
  geom_line(aes(x = Year, y = V, colour = Scenario))



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


g_avt


sim0 <- read_csv(here::here("out", "post_dyage", "Sim_IntvAll.csv"))
tab_avt0 <- sim0 %>% filter(Year >= 2024 & Year <= 2040) %>% 
  select(Year, CumInc, CumMor, Scenario, Key) %>% 
  left_join(
    sim0  %>% 
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

g_avtinc <- tab_avt %>% 
  mutate(
    Scenario = scs[Scenario],
    Scenario = factor(Scenario, scs),
  ) %>% 
  filter(name == "AvtInc") %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
  geom_line(aes(x = Year, y = M, colour = Scenario)) +
  geom_line(data = tab_avt0 %>% filter(name == "AvtInc"), aes(x = Year, y = M, colour = Scenario), linetype = 2) +
  scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
  scale_y_continuous("Averted incidence, %", labels = scales::percent) +
  #facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
  expand_limits(y = 0)

g_avtinc

g_avtmor <- tab_avt %>% 
  mutate(
    Scenario = scs[Scenario],
    Scenario = factor(Scenario, scs)
  ) %>% 
  filter(name == "AvtMor") %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
  geom_line(aes(x = Year, y = M, colour = Scenario)) +
  geom_line(data = tab_avt0 %>% filter(name == "AvtMor"), aes(x = Year, y = M, colour = Scenario), linetype = 2) +
  scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
  scale_y_continuous("Averted mortality, %", labels = scales::percent) +
  #facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
  expand_limits(y = 0)

g_avtmor

ggsave(g_avtinc, filename = here::here("docs", "figs", "g_intv_tx_avtinc.png"), width = 6, height = 4)
ggsave(g_avtmor, filename = here::here("docs", "figs", "g_intv_tx_avtmor.png"), width = 6, height = 4)



sim <- read_csv(here::here("out", folder, "Sim_IntvTx_RelRed.csv"))
sim <- sim %>% filter(Year >= 2024 & Year <= 2040)

sim %>% 
  group_by(Scenario, Year) %>% 
  summarise(V = mean(IncR)) %>% 
  ggplot() +
  geom_line(aes(x = Year, y = V, colour = Scenario))



tab_avt <- sim %>% 
  select(Year, CumInc, CumMor, Scenario, Key) %>% 
  left_join(
    sim  %>% 
      filter(Scenario == "RelRed_1") %>% 
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



g_grad <- tab_avt  %>% 
  filter(Scenario != "RelRed_1") %>% 
  mutate(
    Scenario = gsub("RelRed_", "", Scenario),
    Red = as.numeric(Scenario)
  ) %>% 
  filter(Year == 2040) %>%
  ggplot() +
  geom_line(aes(x = Red, y = M)) + 
  geom_pointinterval(aes(x = Red, y = M, ymax = U, ymin = L)) +
  scale_x_continuous("Percentage reduction in relapse rate, %", labels = scales::percent) + 
  scale_y_continuous("Averted cases, 2024-2040, %", labels = scales::percent) +
  facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
  expand_limits(y = 0)

g_grad


g_gradinc <- tab_avt  %>% 
  filter(Scenario != "RelRed_1") %>% 
  mutate(
    Scenario = gsub("RelRed_", "", Scenario),
    Red = as.numeric(Scenario)
  ) %>% 
  filter(name == "AvtInc") %>% 
  filter(Year == 2040) %>%
  ggplot() +
  geom_line(aes(x = Red, y = M)) + 
  geom_pointinterval(aes(x = Red, y = M, ymax = U, ymin = L)) +
  scale_x_continuous("Percentage reduction in relapse rate, %", labels = scales::percent) + 
  scale_y_continuous("Averted incident cases, 2024-2040, %", labels = scales::percent) +
  expand_limits(y = 0)

g_gradinc



ggsave(g_grad, filename = here::here("docs", "figs", "g_intv_tx_gradpan.png"), width = 6, height = 4)
ggsave(g_gradinc, filename = here::here("docs", "figs", "g_intv_tx_gradincpan.png"), width = 6, height = 4)




