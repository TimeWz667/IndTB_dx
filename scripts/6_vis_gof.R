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


sims %>%
  ggplot(aes(x = Year, y = IncR)) +
  stat_lineribbon(aes(y = IncR), .width = c(.99, .95, .8, .5), color = "#08519C") +
  geom_pointrange(data = tar_inc, aes(x = Year, y = M, ymin = L, ymax = U)) +
  expand_limits(y = 0) +
  scale_y_continuous("Per 100k", labels = scales::number_format(scale = 1e5))
