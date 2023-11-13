library(tidyverse)
library(tidybayes)


theme_set(theme_bw())


qci <- 0.5

## Loading data -----
sim <- read_csv(here::here("out", "dyage", "Sim_IntvAllasassd.csv"))
sim <- sim %>% filter(Year >= 2024 & Year <= 2040)


sim %>% 
  filter(Key == 0) %>% 
  ggplot() +
  geom_line(aes(x = Year, y = IncRecentR, colour = Scenario, group = paste(Scenario, Key)))

avt <- sim %>% 
  select(Year, CumInc, CumMor, Scenario, Key) %>% 
  left_join(
    sim  %>% 
      filter(Scenario == "Baseline") %>% 
      select(Year, CumInc0 = CumInc, CumMor0 = CumMor, Key)
  ) %>% 
  group_by(Key, Scenario) %>% 
  mutate(
    CumInc = CumInc - CumInc[1],
    CumMor = CumMor - CumMor[1],
    CumInc0 = CumInc0 - CumInc0[1],
    CumMor0 = CumMor0 - CumMor0[1]
  ) %>% 
  ungroup() %>% 
  mutate(
    DiffInc = CumInc0 - CumInc,
    DiffMor = CumMor0 - CumMor,
    AvtInc = DiffInc / CumInc0,
    AvtMor = DiffMor / CumMor0,
    AvtInc = ifelse(is.na(AvtInc), 0, AvtInc),
    AvtMor = ifelse(is.na(AvtMor), 0, AvtMor)
  )


avt_rr <- sim %>% 
  select(Year, CumIncRecent, CumIncRemote, Scenario, Key) %>% 
  left_join(
    sim  %>% 
      filter(Scenario == "Baseline") %>% 
      select(Year, CumIncRecent0 = CumIncRecent, CumIncRemote0 = CumIncRemote, Key)
  ) %>% 
  group_by(Key, Scenario) %>% 
  mutate(
    CumIncRecent = CumIncRecent - CumIncRecent[1],
    CumIncRemote = CumIncRemote - CumIncRemote[1],
    CumIncRecent0 = CumIncRecent0 - CumIncRecent0[1],
    CumIncRemote0 = CumIncRemote0 - CumIncRemote0[1]
  ) %>% 
  ungroup() %>% 
  mutate(
    DiffIncRc = CumIncRecent0 - CumIncRecent,
    DiffIncRm = CumIncRemote0 - CumIncRemote,
    AvtIncRc = DiffIncRc / CumIncRecent0,
    AvtIncRm = DiffIncRm / CumIncRemote0,
    AvtIncRc = ifelse(is.na(AvtIncRc), 0, AvtIncRc),
    AvtIncRm = ifelse(is.na(AvtIncRm), 0, AvtIncRm)
  )



tab_epi <- sim %>% 
  select(Year, IncR, MorR, Scenario) %>% 
  pivot_longer(c(IncR, MorR)) %>% 
  group_by(Year, Scenario, name) %>% 
  summarise(
    M = median(value),
    L = quantile(value, (1 - qci) / 2),
    U = quantile(value, 1 - (1 - qci) / 2)
  )


tab_avt <- avt %>% 
  select(Year, AvtInc, AvtMor, Scenario) %>% 
  pivot_longer(c(AvtInc, AvtMor)) %>% 
  group_by(Year, Scenario, name) %>% 
  summarise(
    M = median(value),
    L = quantile(value, (1 - qci) / 2),
    U = quantile(value, 1 - (1 - qci) / 2)
  )


tab_epi_rr <- sim %>% 
  select(Year, IncRcR = IncRecentR, IncRmR = IncRemoteR, Scenario) %>% 
  pivot_longer(c(IncRcR, IncRmR)) %>% 
  group_by(Year, Scenario, name) %>% 
  summarise(
    M = median(value),
    L = quantile(value, (1 - qci) / 2),
    U = quantile(value, 1 - (1 - qci) / 2)
  )

tab_avt_rr <- avt_rr %>% 
  select(Year, AvtIncRc, AvtIncRm, Scenario) %>% 
  pivot_longer(c(AvtIncRc, AvtIncRm,)) %>% 
  group_by(Year, Scenario, name) %>% 
  summarise(
    M = median(value),
    L = quantile(value, (1 - qci) / 2),
    U = quantile(value, 1 - (1 - qci) / 2)
  )



write_csv(tab_epi, here::here("docs", "tabs", "intv_epi.csv"))
write_csv(tab_avt, here::here("docs", "tabs", "intv_avt.csv"))


fn_plot <- function(tab_avt, tab_epi, tab_avt_rr, tab_epi_rr, scs) {
  g_epi <- tab_epi %>%
    filter(Scenario %in% names(scs)) %>% 
    mutate(
      PPM = ifelse(str_detect(Scenario, "PPM"), "With PPM", "Without PPM"),
      Scenario = gsub("PPM", "", Scenario)
    ) %>% 
    mutate(
      Scenario = scs[Scenario],
      Scenario = factor(Scenario, scs)
    ) %>% 
    ggplot() +
    geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
    geom_line(aes(x = Year, y = M, colour = Scenario)) +
    scale_y_continuous("per 100k", labels = scales::number_format(scale = 1e5)) +
    scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
    facet_wrap(.~name, scale = "free_y", 
               labeller = labeller(name = c("IncR" = "Incidence", "MorR" = "Mortality"))) +
    expand_limits(y = 0)
  
  
  g_avt <- tab_avt %>% 
    filter(Scenario %in% names(scs)) %>% 
    mutate(
      PPM = ifelse(str_detect(Scenario, "PPM"), "With PPM", "Without PPM"),
      Scenario = gsub("PPM", "", Scenario)
    ) %>% 
    mutate(
      Scenario = scs[Scenario],
      Scenario = factor(Scenario, scs)
    ) %>% 
    ggplot() +
    geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
    geom_line(aes(x = Year, y = M, colour = Scenario)) +
    scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
    scale_y_continuous("per cent", labels = scales::percent) +
    facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
    expand_limits(y = 0)
  
  
  g_epi_rr <- tab_epi_rr %>%
    filter(Scenario %in% names(scs)) %>% 
    mutate(
      PPM = ifelse(str_detect(Scenario, "PPM"), "With PPM", "Without PPM"),
      Scenario = gsub("PPM", "", Scenario)
    ) %>% 
    mutate(
      Scenario = scs[Scenario],
      Scenario = factor(Scenario, scs)
    ) %>% 
    ggplot() +
    geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
    geom_line(aes(x = Year, y = M, colour = Scenario)) +
    scale_y_continuous("per 100k", labels = scales::number_format(scale = 1e5)) +
    scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
    facet_wrap(.~name, scale = "free_y", 
               labeller = labeller(name = c("IncRcR" = "Incidence, recent", "IncRmR" = "Incidence, remote"))) +
    expand_limits(y = 0)
  
  
  g_avt_rr <- tab_avt_rr %>% 
    filter(Scenario %in% names(scs)) %>% 
    mutate(
      PPM = ifelse(str_detect(Scenario, "PPM"), "With PPM", "Without PPM"),
      Scenario = gsub("PPM", "", Scenario)
    ) %>% 
    mutate(
      Scenario = scs[Scenario],
      Scenario = factor(Scenario, scs)
    ) %>% 
    ggplot() +
    geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
    geom_line(aes(x = Year, y = M, colour = Scenario)) +
    scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
    scale_y_continuous("per cent", labels = scales::percent) +
    facet_wrap(.~name, labeller = labeller(name = c("AvtIncRc" = "Averted recent cases", "AvtIncRm" = "Averted remote cases"))) +
    expand_limits(y = 0)
  
  list(
    g_epi = g_epi,
    g_avt = g_avt,
    g_epi_rr = g_epi_rr,
    g_avt_rr = g_avt_rr
  )
}



gs_dx <- fn_plot(tab_avt, tab_epi, tab_avt_rr, tab_epi_rr, c(
  Baseline = "Baseline",
  Dx_TSwab = "Swab with current tests",
  Dx_POC = "Swab with PoC test",
  Dx_POC_Hi = "Swab with\nhigh-throughput NAAT"
))

gs_dx_itt <- fn_plot(tab_avt, tab_epi, tab_avt_rr, tab_epi_rr, c(
  Baseline = "Baseline",
  Dx_TSwab = "Swab with current tests",
  Dx_POC = "Swab with PoC test",
  Dx_POC_Hi = "Swab with\nhigh-throughput NAAT",
  Dx_TSwab_ITT = "Swab with current tests \nhalf those without intention to test",
  Dx_POC_ITT = "Swab with PoC test \nhalf those without intention to test",
  Dx_POC_Hi_ITT = "Swab with\nhigh-throughput NAAT \nhalf those without intention to test"
))

gs_tx <- fn_plot(tab_avt, tab_epi, tab_avt_rr, tab_epi_rr, c(
  Baseline = "Baseline",
  "Tx_PAN-TB" = "PAN-TB",
  "Tx_LA-INJ" = "LA-INJ"
))

gs_vac <- fn_plot(tab_avt, tab_epi, tab_avt_rr, tab_epi_rr, c(
  Baseline = "Baseline",
  "Vac_BCG" = "BCG revaccination",
  "Vac_M72" = "M72",
  "Vac_BCG-M72" = "Both",
  "Vac_Recurrence" = "Prevent recurrence"
))


gs_mass <- fn_plot(tab_avt, tab_epi, tab_avt_rr, tab_epi_rr, c(
  Baseline = "Baseline",
  Mass_NAAT_20_HRZE = "Mass screening 20% CB-NAAT",
  Mass_NAAT_20_PAN  = "Mass screening 20% CB-NAAT \nwith PAN-TB uptake",
  Mass_Xray_10_HRZE = "Mass screening 10% CXR",
  Mass_Xray_10_PAN = "Mass screening 10% CXR \nwith PAN-TB uptake"
))


gs_combine <- fn_plot(tab_avt, tab_epi, tab_avt_rr, tab_epi_rr, c(
  Baseline = "Baseline",
  "Combine_Lo" = "Combination (conservative)",
  "Combine_Hi" = "Combination (ambitious)"
))

# 
# scs <- c(
#   Baseline = "Baseline",
#   "CombinePPM_Lo" = "Combination A",
#   "CombinePPM_Hi" = "Combination B",
#   "Combine_Lo" = "Combination A",
#   "Combine_Hi" = "Combination B"
# )
# 
# g_combine <- tab_avt %>% 
#   ungroup() %>% 
#   filter(Scenario %in% names(scs)) %>% 
#   mutate(
#     PPM = ifelse(str_detect(Scenario, "PPM"), "With PPM", "Without PPM"),
#     Scenario = gsub("PPM", "", Scenario)
#   ) %>% 
#   mutate(
#     Scenario = scs[Scenario],
#     Scenario = factor(Scenario, unique(scs))
#   ) %>% 
#   ggplot() +
#   #geom_ribbon(aes(x = Year, ymin=L, ymax=U, fill=Scenario), alpha = 0.1) +
#   geom_line(aes(x = Year, y = M, colour = Scenario, linetype = PPM)) +
#   scale_x_continuous("Year", breaks = c(2024, seq(2030, 2040, 5))) + 
#   scale_y_continuous("per cent", labels = scales::percent) +
#   facet_wrap(.~name, labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths"))) +
#   expand_limits(y = 0)



ggsave(gs_dx$g_epi, filename = here::here("docs", "figs", "g_intv_dx_epi.png"), width = 12, height = 5)
ggsave(gs_dx$g_avt, filename = here::here("docs", "figs", "g_intv_dx_avt.png"), width = 12, height = 5)

ggsave(gs_dx_itt$g_epi, filename = here::here("docs", "figs", "g_intv_dx_itt_epi.png"), width = 12, height = 5)
ggsave(gs_dx_itt$g_avt, filename = here::here("docs", "figs", "g_intv_dx_itt_avt.png"), width = 12, height = 5)


ggsave(gs_tx$g_epi, filename = here::here("docs", "figs", "g_intv_tx_epi.png"), width = 12, height = 5)
ggsave(gs_tx$g_avt, filename = here::here("docs", "figs", "g_intv_tx_avt.png"), width = 12, height = 5)

ggsave(gs_vac$g_epi, filename = here::here("docs", "figs", "g_intv_vac_epi.png"), width = 12, height = 5)
ggsave(gs_vac$g_avt, filename = here::here("docs", "figs", "g_intv_vac_avt.png"), width = 12, height = 5)

ggsave(gs_mass$g_epi, filename = here::here("docs", "figs", "g_intv_acf_epi.png"), width = 12, height = 5)
ggsave(gs_mass$g_avt, filename = here::here("docs", "figs", "g_intv_acf_avt.png"), width = 12, height = 5)

ggsave(gs_combine$g_epi, filename = here::here("docs", "figs", "g_intv_combine_epi.png"), width = 12, height = 5)
ggsave(gs_combine$g_avt, filename = here::here("docs", "figs", "g_intv_combine_avt.png"), width = 12, height = 5)



scs <- c(
  "Tx_PAN-TB" = "PAN-TB",
  "Tx_LA-INJ" = "LA-INJ",
  Dx_TSwab = "Swab with current tests",
  Dx_POC = "Swab with PoC test",
  Dx_POC_Hi = "Swab with\nhigh-throughput CB-NAAT",
  Mass_NAAT_20_HRZE = "Mass screening 20% CB-NAAT",
  Mass_Xray_10_HRZE = "Mass screening 10% CXR",
  "Vac_BCG" = "BCG revaccination",
  "Vac_M72" = "M72",
  "Vac_BCG-M72" = "BCG revac + M72",
  "Vac_Recurrence" = "Prevent recurrence"
)


g_intv_summary <- tab_avt %>% 
  ungroup() %>% 
  filter(Year == 2040) %>% 
  filter(!startsWith(Scenario, "Combine")) %>% 
  mutate(
    PPM = ifelse(str_detect(Scenario, "PPM"), "With PPM", "Without PPM"),
    Scenario = gsub("PPM", "", Scenario)
  ) %>% 
  filter(Scenario %in% names(scs)) %>% 
  mutate(
    Scenario = scs[Scenario],
    Scenario = factor(Scenario, scs)
  ) %>% 
  ggplot() +
  geom_bar(aes(x = Scenario, y = M, fill = PPM), stat = "identity", position = "dodge") +
  scale_y_continuous("percent cases averted", labels = scales::percent) +
  scale_fill_discrete("") +
  theme(axis.text.x.bottom = element_text(hjust = 0, angle = -30)) +
  facet_grid(name~.,labeller = labeller(name = c("AvtInc" = "Averted cases", "AvtMor" = "Averted deaths")))

ggsave(g_intv_summary, filename = here::here("docs", "figs", "g_intv_summary.png"), width = 12, height = 5)

