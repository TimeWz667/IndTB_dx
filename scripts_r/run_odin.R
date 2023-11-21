library(odin)
library(tidyverse)

theme_set(theme_bw())

pars_demo <- jsonlite::read_json(here::here("pars", "ind_who_90to45.json"), simplifyVector = T)

dim(pars_demo$RateDeath)


model <- odin::odin("model/tb.R", target = "c")


pars <- with(pars_demo,{
  list(
    tt_demo = Year,
    n0 = N,
    n_agp = dim(RateDeath)[2],
    r_death = RateDeath,
    r_aging = RateAgeing,
    r_birth = RateBirth
  )
})

pars$Y0 <- array(0, c(16, 3, 7))
pars$Y0[1, 1, ] <- pars$n0[1, ]


cm <- model$new(user = pars)

ys <- cm$run(seq(1990, 2041, 0.5))




plot(ys[, "N"], type = 'l')
plot(pars$tt_demo, rowSums(pars$n0), cex = 1.5)
lines(ys[, "t"], ys[, "N"])


ns <- as_tibble(as.data.frame(ys)) %>% 
  select(t, starts_with("Y[")) %>% 
  pivot_longer(-t, values_to = "N") %>%
  extract(name, c("TB", "DR", "Agp"), "Y\\[(\\d+),(\\d+),(\\d+)\\]", convert = T) %>% 
  mutate(
    Agp = pars_demo$dimnames$Age[Agp],
    Agp = factor(Agp, pars_demo$dimnames$Age)
  )


dat_ns <- data.frame(t = pars_demo$Year - 0.5, pars_demo$N) %>% 
  pivot_longer(-t, values_to = "N_Data") %>%
  extract(name, c("Agp"), "X(\\d+)", convert = T) %>% 
  mutate(
    Agp = pars_demo$dimnames$Age[Agp],
    Agp = factor(Agp, pars_demo$dimnames$Age)
  )



g_pop <- ns %>% 
  filter(t >= 2000) %>% 
  group_by(t) %>% 
  summarise(N = sum(N))%>% 
  left_join(dat_ns %>% 
              mutate(t = t + 0.5) %>% 
              group_by(t) %>% 
              summarise(N_Data = sum(N_Data))) %>% 
  ggplot() +
  geom_line(aes(x = t, y = N)) +
  geom_point(aes(x = t, y = N_Data, colour = "Data")) +
  scale_y_continuous("Population, Million", labels = scales::number_format(scale = 1e-6)) +  
  scale_x_continuous("Year") +
  scale_color_discrete("") +
  expand_limits(y = 0) +
  theme(legend.position = c(1, 0), legend.justification = c(1.05, -0.05))
  



g_agp <- ns %>% 
  filter(t %in% seq(2000, 2040, 10)) %>% 
  left_join(dat_ns) %>% 
  ggplot() +
  geom_bar(aes(x = N, y = Agp), alpha = 0.3, stat = "identity") +
  geom_point(aes(x = N_Data, y = Agp, colour = "Data")) +
  scale_x_continuous("Population, Million", labels = scales::number_format(scale = 1e-6)) +  
  scale_color_discrete("") +
  facet_wrap(.~t)

g_pop
g_agp


ggsave(g_pop, filename = here::here("docs", "figs", "g_baseline_pop.png"), width = 7, height = 5)
ggsave(g_agp, filename = here::here("docs", "figs", "g_baseline_popy.png"), width = 7, height = 5)

