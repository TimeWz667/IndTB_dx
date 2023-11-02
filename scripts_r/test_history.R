library(tidyverse)


theme_set(theme_bw())


ms <- read_csv(here::here("out", "sim_ppv_treated.csv"))

g <- ms %>% 
  ggplot() + 
  geom_line(aes(x = risk, y = pi, colour = "True History")) +
  geom_line(aes(x = risk, y = ppv, colour = "True Cases")) +
  scale_y_continuous("Proportional true TB, %", labels = scales::percent) +
  scale_x_continuous("Life-time risk of TB infection, %", labels = scales::percent) +
  scale_color_discrete("") + 
  guides(colour = guide_legend(reverse=T)) + 
  expand_limits(y = 0:1)

g


ggsave(g, filename="out/test.png", width = 5, height = 4)
