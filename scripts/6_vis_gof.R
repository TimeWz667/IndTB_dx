library(tidyverse)
library(tidybayes)

theme_set(theme_bw())



sims <- read_csv(here::here("out", "post_dy", "Sim.csv"))

post <- 
  read_csv(here::here("out", "post_dy", "post.csv"))

head(sims)

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

