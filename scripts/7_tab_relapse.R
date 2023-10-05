library(tidyverse)



post <- read_csv(here::here("out", "post_dyage", "Post.csv"))


summ <- function(df) {
  df %>% 
    mutate(
      rr_rel_pub = rr_relapse_pub / (0.3 + 0.7 * rr_relapse_pub),
      rr_rel_pri = 1 / (0.3 + 0.7 * rr_relapse_pub),
      r_rel_pub = r_relapse_te * k_relapse_adj * rr_rel_pub,
      r_rel_pri = r_relapse_te * k_relapse_adj * rr_rel_pri,
      p_pub = pexp(1, r_rel_pub),
      p_pri = pexp(1, r_rel_pri)
    ) %>% 
    summarise(across(r_rel_pub:p_pri, list(
      m = median,
      l = function(x) quantile(x, 0.025),
      u = function(x) quantile(x, 0.975)
    ))) %>% 
    pivot_longer(everything()) %>% 
    extract(name, c("Index", "name"), "(\\S+)_(m|l|u)") %>% 
    pivot_wider()
}


post %>% 
  select(r_relapse_te, rr_relapse_pub) %>% 
  mutate(
    k_relapse_adj = 2 / 1.2
  ) %>% 
  summ()


post %>% 
  select(r_relapse_te, rr_relapse_pub) %>% 
  mutate(
    k_relapse_adj = 3.2
  ) %>% 
  summ()


post %>% 
  select(r_relapse_te, rr_relapse_pub) %>% 
  mutate(
    k_relapse_adj = 4.4,
    rr_relapse_pub = 1
  ) %>% 
  summ()

