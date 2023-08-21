library(tidyverse)

theme_set(theme_bw())


targets <- read_csv(here::here("data", "targets.csv")) %>% 
  filter(Year >= 2016) %>%
  filter(Year <= 2020) %>% 
  filter(Index %in% c("Inc", "Mor")) %>% 
  mutate(
    Index = ifelse(Index == "Inc", "IncR", "MorR")
  )
mss0 <- read_csv(here::here("out", "intv", "mss_pre.csv"))
mss1 <- read_csv(here::here("out", "intv", "mss_intv.csv"))


scs <- c(
  Baseline = "Baseline",
  Without_LTBI = "Impacts on Active TB early detection",
  With_LTBI = "Considering impacts on LTBIs"
)



#### Trends on epidemiology -----

epi0 <- mss0 %>% 
  select(Year, IncR, MorR) %>% 
  pivot_longer(-Year, names_to = "Index") %>% 
  group_by(Year, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  )


epi1 <- mss1 %>% 
  select(Year, Scenario, IncR, MorR) %>% 
  pivot_longer(-c(Year, Scenario), names_to = "Index") %>% 
  group_by(Year, Index, Scenario) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  )


epi <- bind_rows(
  epi0 %>% 
    filter(Year >= 2016) %>% 
    mutate(Scenario = "Baseline") %>% 
    ungroup(),
  epi1  %>% 
    filter(Year > 2022)%>% 
    filter(Year <= 2030)%>% 
    ungroup() 
) %>% 
  ungroup() %>% 
  mutate(
    Scenario = factor(Scenario, names(scs))
  )


g_epi <- epi %>% 
  ggplot(aes(x = Year)) + 
  geom_ribbon(aes(ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(y = M, colour = Scenario)) +
  geom_pointrange(data = targets, aes(y = M, ymin = L, ymax = U)) + 
  scale_fill_discrete(labels = scs) +
  scale_colour_discrete(labels = scs) +
  scale_x_continuous(breaks = c(2016, 2020, 2023, 2025, 2030)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(.~Index, scale = "free_y", labeller = labeller(Index=c(IncR="Incidence", MorR = "Mortality"))) +
  expand_limits(y = 0) +
  geom_vline(xintercept = 2023, linetype = 2) +
  theme(legend.position = "bottom")

g_epi


#### Incremental aversion -----

avt <- local({
  avt <- mss1 %>% 
    filter(Year >= 2023) %>% 
    filter(Year <= 2030) %>% 
    select(Year, Scenario, Key, CumInc, CumMor) %>% 
    pivot_longer(-c(Year, Scenario, Key), names_to = "Index") %>% 
    group_by(Scenario, Key, Index) %>% 
    arrange(Year) %>% 
    mutate(value = value - value[1]) %>% 
    ungroup()
  
  avt <- avt %>%
    left_join(avt %>%  filter(Scenario == "Baseline") %>% select(-Scenario, value0 = value)) %>% 
    mutate(value = - (value / value0 - 1), value = ifelse(is.na(value), 0, value)) %>% 
    group_by(Year, Index, Scenario) %>% 
    summarise(
      M = median(value),
      L = quantile(value, 0.025),
      U = quantile(value, 0.975)
    ) %>% 
    mutate(
      Scenario = factor(Scenario, names(scs))
    )
  
  avt
})


g_avt <- avt %>% 
  ggplot(aes(x = Year)) + 
  geom_ribbon(aes(ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(y = M, colour = Scenario)) +
  scale_fill_discrete(labels = scs) +
  scale_colour_discrete(labels = scs) +
  scale_x_continuous(breaks = c(2016, 2020, 2023, 2025, 2030)) +
  scale_y_continuous("percentage averted", labels = scales::percent) + 
  facet_wrap(.~Index, scale = "free_y", labeller = labeller(Index=c(CumInc="Incidence", CumMor = "Mortality"))) +
  expand_limits(y = 0) +
  geom_vline(xintercept = 2023, linetype = 2) +
  theme(legend.position = "bottom")

g_avt



#### Unnecessary treatment -----

fp <- mss1 %>% 
  filter(Year >= 2023) %>% 
  filter(Year <= 2030) %>% 
  select(Year, Scenario, Key, CumYieldLTBI, CumYieldOT) %>% 
  pivot_longer(-c(Year, Scenario, Key), names_to = "Index") %>% 
  group_by(Scenario, Key, Index) %>% 
  arrange(Year) %>% 
  mutate(value = value - value[1]) %>% 
  group_by(Year, Index, Scenario) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  ) %>% 
  mutate(
    Scenario = factor(Scenario, names(scs))
  )


g_fp <- fp %>% 
  ggplot(aes(x = Year)) + 
  geom_ribbon(aes(ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(y = M, colour = Scenario)) +
  scale_fill_discrete(labels = scs) +
  scale_colour_discrete(labels = scs) +
  scale_x_continuous(breaks = c(2016, 2020, 2023, 2025, 2030)) +
  scale_y_continuous("Yields via false-positivity, million", labels = scales::number_format(scale = 1e-6)) + 
  facet_wrap(.~Index, scale = "free_y", labeller = labeller(Index=c(CumYieldLTBI="Latent TB", CumYieldOT="Uninfected and successfully treated"))) +
  expand_limits(y = 0) +
  geom_vline(xintercept = 2023, linetype = 2) +
  theme(legend.position = "bottom")

g_fp


#### LTBI -----

ltbi0 <- mss0 %>% 
  mutate(Recent = IncRecentR / IncR, Retreated = IncRetreatR / IncR) %>% 
  select(Year, LTBI, Recent, Retreated) %>% 
  pivot_longer(-Year, names_to = "Index") %>% 
  group_by(Year, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  )


ltbi1 <- mss1 %>% 
  mutate(Recent = IncRecentR / IncR, Retreated = IncRetreatR / IncR) %>% 
  select(Year, Scenario, LTBI, Recent, Retreated) %>% 
  pivot_longer(-c(Year, Scenario), names_to = "Index") %>% 
  group_by(Year, Index, Scenario) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  )


ltbi <- bind_rows(
  ltbi0 %>% 
    filter(Year >= 2016) %>% 
    mutate(Scenario = "Baseline") %>% 
    ungroup(),
  # ltbi1  %>% 
  #   filter(Year > 2022)%>% 
  #   filter(Year <= 2030)%>% 
  #   ungroup() 
) %>% 
  ungroup() %>% 
  mutate(
    Scenario = factor(Scenario, names(scs))
  )


g_ltbi <- ltbi %>% 
  ggplot(aes(x = Year)) + 
  geom_ribbon(aes(ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(y = M, colour = Scenario)) +
  scale_fill_discrete(labels = scs) +
  scale_colour_discrete(labels = scs) +
  scale_x_continuous(breaks = c(2016, 2020, 2023, 2025, 2030)) +
  scale_y_continuous("percent", labels = scales::percent) + 
  facet_wrap(.~Index, scale = "free_y", labeller = labeller(Index=c(LTBI="Latent TB infection", 
                                                                    Recent = "Incident TB \n from recent infection",
                                                                    Retreated = "Incident TB \nwith treatment history"))) +
  expand_limits(y = 0:1) +
  #geom_vline(xintercept = 2023, linetype = 2) +
  theme(legend.position = "bottom")

g_ltbi



#### Output -----


ggsave(g_epi, filename=here::here("docs", "figs", "g_bs_epi.png"), width = 8, height = 5)
ggsave(g_avt, filename=here::here("docs", "figs", "g_bs_avt.png"), width = 8, height = 5)
ggsave(g_fp, filename=here::here("docs", "figs", "g_bs_fp.png"), width = 8, height = 5)
ggsave(g_ltbi, filename=here::here("docs", "figs", "g_bs_ltbi.png"), width = 8, height = 5)





mss1 %>% 
  select(Year, Key, Scenario, NotiR, IncR, YieldATBR) %>% 
  group_by(Scenario) %>% 
  mutate(YieldATBR / NotiR)

mss1 %>% 
  select(Year, Key, Scenario, NotiR, IncR, YieldATBR) %>% 
  left_join(mss1 %>% filter(Scenario == "Baseline") %>% select(Year, Key, NotiR0 = NotiR)) %>% 
  mutate(PrIncre = YieldATBR / NotiR0, Add = (YieldATBR + NotiR) / NotiR0) %>% 
  group_by(Year, Scenario) %>% 
  summarise(
    Pi_M = median(PrIncre), 
    Pi_L = quantile(PrIncre, 0.025),
    Pi_U = quantile(PrIncre, 0.975),
    Add_M = median(Add), 
    Add_L = quantile(Add, 0.025),
    Add_U = quantile(Add, 0.975)
  ) %>%  
  filter(Year == 2030)


mss1 %>% 
  mutate(
    Dur = PrevUt / IncR
  ) %>% 
  select(Year, Key, Scenario, Dur) %>% 
  group_by(Year, Scenario) %>% 
  summarise(
    Dur_M = median(Dur), 
    Dur_L = quantile(Dur, 0.025),
    Dur_U = quantile(Dur, 0.975)
  ) %>% 
  arrange(Scenario) %>% 
  data.frame()
