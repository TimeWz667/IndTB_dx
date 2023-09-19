library(tidyverse)
library(tidybayes)

theme_set(theme_bw())



sims <- read_csv(here::here("out", "post_dy", "Sim.csv"))

targets <- read_csv(here::here("data", "targets.csv"))


tar_inc <- targets %>% 
  filter(Year >= 2015) %>% 
  filter(Tag == "All" & Index == "IncR")


g_inc <- sims %>%
  ggplot(aes(x = Year, y = IncR)) +
  stat_lineribbon(aes(y = IncR), .width = c(.99, .95, .8, .5), color = "#08519C") +
  geom_pointrange(data = tar_inc, aes(x = Year, y = M, ymin = L, ymax = U)) +
  expand_limits(y = 0) +
  scale_y_continuous("Incidence per 100k", labels = scales::number_format(scale = 1e5)) +
  scale_fill_brewer() 

ggsave(g_inc, filename = here::here("docs", "figs", "g_baseline_inc.png"), width = 7, height = 5)



g_treated <- sims %>% 
  mutate(PrPub = IncTreatedPubR / IncTreatedR) %>%
  ggplot(aes(x = Year, y = PrPub)) +
  stat_lineribbon(aes(y = PrPub), .width = c(.99, .95, .8, .5), color = "#08519C") +
  geom_hline(yintercept = 0.81) +
  expand_limits(y = 0:1) +
  scale_y_continuous("Previously treated in public sectors / Cases with treatment history, Percent", labels = scales::percent) +
  scale_fill_brewer() 

ggsave(g_treated, filename = here::here("docs", "figs", "g_baseline_treated.png"), width = 7, height = 5)




g_inc_sub <- sims %>% 
  select(Year, starts_with("Inc")) %>% 
  pivot_longer(-Year) %>% 
  group_by(Year, name) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = value, colour = name)) +
  expand_limits(y = 0) +
  scale_y_continuous("Incidence per 100k", labels = scales::number_format(scale = 1e5)) +
  scale_colour_discrete("Source of cases", labels = c(
    "IncR" = "All",
    "IncRecentR" = "Infected < 2yr",
    "IncRemoteR" = "Infected > 2yr",
    "IncTreatedPubR" = "Treated in public hosp.",
    "IncTreatedR" = "Treated"
  )) 


ggsave(g_inc_sub, filename = here::here("docs", "figs", "g_baseline_inc_sub.png"), width = 7, height = 5)



g_ltbi <- sims %>% 
  ggplot(aes(x = Year)) +
  stat_lineribbon(aes(y = LTBI), .width = c(.99, .95, .8, .5), color = "#08519C") +
  expand_limits(y = 0) +
  scale_y_continuous("LTBI, %", labels = scales::percent) +
  scale_fill_brewer() 

ggsave(g_ltbi, filename = here::here("docs", "figs", "g_baseline_ltbi.png"), width = 7, height = 5)



## Age group model -----
sims <- read_csv(here::here("out", "post_dyage", "Sim.csv"))

targets <- read_csv(here::here("data", "targets.csv"))


tar_inc <- targets %>% 
  filter(Year >= 2015) %>% 
  filter(Tag == "All" & Index == "IncR")

tar_inca <- targets %>% 
  filter(Year >= 2015) %>% 
  filter(Tag != "All" & Index == "IncR") %>% 
  filter(Tag != "0-14") %>% 
  rename(AgeGrp = Tag)


g_inc <- sims %>%
  ggplot(aes(x = Year, y = IncR)) +
  stat_lineribbon(aes(y = IncR), .width = c(.99, .95, .8, .5), color = "#08519C") +
  geom_pointrange(data = tar_inc, aes(x = Year, y = M, ymin = L, ymax = U)) +
  expand_limits(y = 0) +
  scale_y_continuous("Incidence per 100k", labels = scales::number_format(scale = 1e5)) +
  scale_fill_brewer() 

g_inc

ggsave(g_inc, filename = here::here("docs", "figs", "g_agp_inc.png"), width = 7, height = 5)



g_treated <- sims %>% 
  mutate(PrPub = IncTreatedPubR / IncTreatedR) %>%
  ggplot(aes(x = Year, y = PrPub)) +
  stat_lineribbon(aes(y = PrPub), .width = c(.99, .95, .8, .5), color = "#08519C") +
  geom_hline(yintercept = 0.81) +
  expand_limits(y = 0:1) +
  scale_y_continuous("Previously treated in public sectors / Cases with treatment history, Percent", labels = scales::percent) +
  scale_fill_brewer() 

g_treated

ggsave(g_treated, filename = here::here("docs", "figs", "g_agp_treated.png"), width = 7, height = 5)




g_inc_sub <- sims %>% 
  select(Year, starts_with("Inc")) %>%
  select(-starts_with("IncR_")) %>% 
  pivot_longer(-Year) %>% 
  group_by(Year, name) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = value, colour = name)) +
  expand_limits(y = 0) +
  scale_y_continuous("Incidence per 100k", labels = scales::number_format(scale = 1e5)) +
  scale_colour_discrete("Source of cases", labels = c(
    "IncR" = "All",
    "IncRecentR" = "Infected < 2yr",
    "IncRemoteR" = "Infected > 2yr",
    "IncTreatedPubR" = "Treated in public hosp.",
    "IncTreatedR" = "Treated"
  )) 

g_inc_sub

ggsave(g_inc_sub, filename = here::here("docs", "figs", "g_agp_inc_sub.png"), width = 7, height = 5)


g_ltbi <- sims %>% 
  ggplot(aes(x = Year)) +
  stat_lineribbon(aes(y = LTBI), .width = c(.99, .95, .8, .5), color = "#08519C") +
  expand_limits(y = 0) +
  scale_y_continuous("LTBI, %", labels = scales::percent) +
  scale_fill_brewer() 

g_ltbi

ggsave(g_ltbi, filename = here::here("docs", "figs", "g_agp_ltbi.png"), width = 7, height = 5)




g_inc_age <- sims %>% 
  select(Year, starts_with("IncR_")) %>%
  pivot_longer(-Year) %>% 
  group_by(Year, name) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = value, colour = name)) +
  expand_limits(y = 0) +
  scale_y_continuous("Incidence per 100k", labels = scales::number_format(scale = 1e5)) +
  scale_colour_discrete("Source of cases", labels = c(
    "IncR" = "All",
    "IncRecentR" = "Infected < 2yr",
    "IncRemoteR" = "Infected > 2yr",
    "IncTreatedPubR" = "Treated in public hosp.",
    "IncTreatedR" = "Treated"
  )) 


sims %>% 
  select(Year, starts_with("IncR_"), starts_with("LTBI_")) %>%
  pivot_longer(-Year) %>% 
  extract(name, c("Index", "AgeGrp"), "(IncR|LTBI)_(\\S+)") %>% 
  filter(!(AgeGrp %in% c("0-4", "5-14"))) %>% 
  filter(Year == 2021) %>% 
  # group_by(AgeGrp, Index) %>% 
  # summarise(
  #   m = median(value),
  #   l = quantile(value, 0.25),
  #   u = quantile(value, 0.75)
  # ) %>% 
  ggplot(aes(y = AgeGrp)) +
  stat_interval(aes(x = value), .width = c(.5, .8, .95)) +
  # geom_histogram(aes(x = m), alpha = 0.3, stat = "identity") +
  # geom_linerange(aes(xmin = l, xmax = u)) + 
  geom_pointinterval(data = tar_inca, aes(x = M, xmin = L, xmax = U), position = position_nudge(y = -0.2)) + 
  scale_x_continuous("percent", labels = scales::percent) + 
  facet_grid(.~Index, scales = "free_x") +
  scale_colour_brewer() +
  expand_limits(x = 0)






