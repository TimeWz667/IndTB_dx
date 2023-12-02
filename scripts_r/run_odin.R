library(odin)
library(tidyverse)

theme_set(theme_bw())

pars_demo <- jsonlite::read_json(here::here("pars", "ind_who_90to45.json"), simplifyVector = T)
pars_cas <- jsonlite::read_json(here::here("pars", "pars_cas_cdx.json"), simplifyVector = T, simplifyDataFrame = F)


pars <- with(pars_demo,{
  list(
    tt_demo = c(0, Year),
    n0 = rbind(N[1, ], N),
    n_agp = dim(RateDeath)[2],
    r_death = rbind(RateDeath[1, ], RateDeath),
    r_aging = rbind(RateAgeing[1, ], RateAgeing),
    r_birth = c(RateBirth[1], RateBirth)
  )
})

pars$Y0 <- array(0, c(16, 3, 7))
pars$Y0[1, 1, ] <- pars$n0[1, ] * 0.95
pars$Y0[2, 1, ] <- pars$n0[1, ] * 0.05

exo <- list(
  r_sc = 0.2,
  rr_die_asym = 1,
  r_die_sym = 0.12
)

p1 <- with(c(exo, pars_cas$Prev, pars_cas$Particles[[1]]), {
  r_die_asym <- r_die_sym * rr_die_asym

  mu_c <- r_die_sym + r_sc
  mu_s <- r_die_sym + r_sc
  mu_a <- r_die_asym + r_sc
  
  p_det <- sum(p_ent * p_itt * p_dx * p_txi)
  
  r_det <- txi / PrevExCS
  r_csi <- (r_det + mu_c) * PrevExCS / PrevSym
  r_onset <- (r_csi + mu_s) * PrevSym / PrevAsym
  inc <- (r_onset + mu_a) * PrevAsym
  
  det0 <- PrevSym * r_csi * p_det
  r_recsi <- (txi - det0) / (PrevExCS * p_det)
  
  list(
    beta = 20,
    irr_25=1,
    irr_35=1,
    irr_45=1,
    irr_55=1,
    irr_65=1,
    r_sc = r_sc,
    rr_die_asym = rr_die_asym,
    r_die_sym = r_die_sym,
    r_onset = r_onset,
    r_csi0 = r_csi,
    r_recsi0 = r_recsi,
    p_ent = p_ent,
    p_itt = p_itt,
    p_dx = p_dx,
    p_txi = p_txi,
    p_pri_on_pub = p_pri_on_pub,
    tx_dur = c(Pu=dur_pub, Pe=dur_pri, Pr=dur_pri),
    p_cure = rbind(
      DS = c(Pu=0.9, Pe=0.85, Pr=0.85),
      DR = c(0.29, 0.29, 0.29),
      NR = c(0.29, 0.29, 0.29)
    )
  )
})


pars_intv <- with(p1, {
  list(
    tt_intv_ms = c(0, 3000),
    uptake_ms = c(0, 0),
    tt_intv_dx = c(0, 3000),
    uptake_dx = c(0, 0),
    tt_intv_tx = c(0, 3000),
    uptake_tx = c(0, 0),
    tt_intv_vac = c(0, 3000),
    uptake_vac = c(0, 0),
    p_txi_new = p_txi,
    p_cure_new = p_cure
  )
})
  


model <- odin::odin("model/tb.R", target = "c")


cm <- model$new(user = c(pars, p1, pars_intv))
ys <- cm$run(seq(1800, 2041, 0.25))
ys <- ys[ys[, "t"] >= 1980, ]

y1 <- array(ys[nrow(ys), startsWith(colnames(ys), "Y[")], c(16, 3, 7))
y1



ys[, "Net_Care"]
ys[, "Net_TB"]
ys[, "Net_Pop"]

plot(ys[, "t"], ys[, "Inc"], type = 'l')

# Plot for demography
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

