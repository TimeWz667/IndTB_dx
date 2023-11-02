library(odin)
library(tidyverse)


pars_demo <- jsonlite::read_json(here::here("pars", "ind_who_70to35.json"), simplifyVector = T)

dim(pars_demo$RateDeath)


model <- odin::odin("R/tb.R", target = "c")


pars <- with(pars_demo,{
  list(
    tt_demo = Year,
    n0 = N,
    r_death = RateDeath,
    r_aging = RateAgeing,
    r_birth = RateBirth
  )
})

pars$Y0 <- array(0, c(16, 3, 8))
pars$Y0[1, 1, ] <- pars$n0[1, ]


cm <- model$new(user = pars)

ys <- cm$run(seq(1970, 2035, 0.5))




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



ns %>% 
  group_by(t) %>% 
  summarise(N = sum(N))%>% 
  left_join(dat_ns %>% 
              mutate(t = t + 0.5) %>% 
              group_by(t) %>% 
              summarise(N_Data = sum(N_Data))) %>% 
  ggplot() +
  geom_line(aes(x = t, y = N)) +
  geom_point(aes(x = t, y = N_Data)) +
  scale_y_continuous("Population, M", labels = scales::number_format(scale = 1e-6)) +  
  scale_x_continuous("Year") +
  expand_limits(y = 0)
  



ns %>% 
  filter(t %in% seq(1970, 2040, 10)) %>% 
  left_join(dat_ns) %>% 
  ggplot() +
  geom_bar(aes(x = N, y = Agp), alpha = 0.3, stat = "identity") +
  geom_point(aes(x = N_Data, y = Agp)) +
  scale_x_continuous("Population, M", labels = scales::number_format(scale = 1e-6)) +  
  facet_wrap(.~t)
  


