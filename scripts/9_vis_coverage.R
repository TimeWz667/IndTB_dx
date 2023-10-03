library(tidyverse)

theme_set(theme_bw())


sim <- read_csv(here::here("out", "coverage_test.csv"))[, -1]

sim <- sim %>% 
  group_by(Coverage, Switch, Cohort, Interval) %>% 
  mutate(
    Inc = c(0, diff(Inc)),
    Mor = c(0, diff(Mor)),
    Covered = c(0, diff(Covered)),
    Eli = c(0, diff(Eli)),
    Yield = c(0, diff(Yield)),
    PYield = c(0, Yield[-1] / Prev[-1])
  ) %>% 
  mutate(Scenario = case_when(
    !is.na(Interval) ~ paste0("Interval: ", Interval),
    !Cohort ~ "r = k",
    T ~ "r = - log(1 - k)"
  )) %>% 
  filter(Switch != "TRUE" | is.na(Switch)) %>% 
  filter(Time > 0)



sim %>% 
  ungroup() %>% 
  mutate(Scenario = paste(Switch, Cohort, Interval)) %>% 
  filter(!is.na(Interval)) %>% 
  ggplot() +
  geom_line(aes(x = Time, y = Inc, colour = as.character(Interval))) +
  scale_y_continuous("Incidence, per 100k", labels = scales::number_format(scale=1e5)) +
  facet_grid(.~Coverage) +
  expand_limits(y = 0)


sim %>% 
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x = Time, y = Inc, colour = Scenario)) +
  scale_y_continuous("Incidence, per 100k", labels = scales::number_format(scale=1e5)) +
  facet_grid(.~Coverage) +
  expand_limits(y = 0)

  
  

sim %>% 
  ungroup() %>% 
  mutate(Scenario = paste(Switch, Cohort, Interval)) %>% 
  #filter(!is.na(Interval)) %>% 
  ggplot() +
  geom_line(aes(x = Time, y = Inc, colour = Scenario)) +
  scale_y_continuous("Incidence, per 100k", labels = scales::number_format(scale=1e5)) +
  facet_grid(.~Coverage) +
  expand_limits(y = 0)


sim %>% 
  #filter(!is.na(Interval)) %>% 
  ggplot() +
  geom_line(aes(x = Time, y = Mor, colour = Scenario)) +
  scale_y_continuous("Mortality, per 100k", labels = scales::number_format(scale=1e5)) +
  facet_grid(.~Coverage) +
  expand_limits(y = 0)



sim %>% 
  ggplot() +
  geom_line(aes(x = Time, y = Yield, colour = Scenario)) +
  scale_y_continuous("Yield", labels = scales::percent) +
  facet_grid(.~Coverage) +
  expand_limits(y = 0)


sim %>% 
  ggplot() +
  geom_line(aes(x = Time, y = PYield, colour = Scenario)) +
  scale_y_continuous("PYield", labels = scales::percent) +
  facet_grid(.~Coverage) +
  expand_limits(y = 0)

