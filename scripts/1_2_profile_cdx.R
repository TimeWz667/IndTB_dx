library(tidyverse)
library(tidybayes)
library(sensitivity)

theme_set(theme_bw())



post <- read_csv(here::here("pars", "pars_cdx.csv"))
post_alt <- read_csv(here::here("pars", "pars_cdx_alt.csv"))


lab.pars <- c(
  "p_tb" = "pr(TB)",
  "p_scanty" = "Sample scanty",
  "sens_s" = "Specificity, SSM",
  "sens_x" = "Specificity, Xpert",
  "sens_x_sn" = "Specificity, Xpert, SSM-",
  "spec_s" = "Specificity, SSM",
  "spec_x" = "Specificity, Xpert",
  "spec_x_sn" = "Specificity, Xpert, SSM-"
)


spec <- pcc(
  X = post %>% select(p_tb, p_scanty, starts_with(c("sens_", "spec_")), -sens_cdx, -spec_cdx, -sens_cdx_bn, -spec_cdx_bn),
  y = post$spec_cdx,
  rank = T
)

g_spec <- tibble(var = rownames(spec$PRCC), est = spec$PRCC$original) %>% 
  arrange(abs(est)) %>% 
  mutate(
    var = factor(var, var)
  ) %>% 
  ggplot() +
  geom_bar(aes(y = var, x = est), stat = "identity") +
  theme(axis.text.x = element_text(angle = -60, vjust = 1, hjust = 0)) +
  scale_x_continuous("PRCC to CDx specificity") +
  scale_y_discrete("Variable", label = lab.pars)


sens <- pcc(
  X = post %>% select(p_tb, p_scanty, starts_with(c("sens_", "spec_")), -sens_cdx, -spec_cdx, -sens_cdx_bn, -spec_cdx_bn),
  y = post$sens_cdx,
  rank = T
)

g_sens <- tibble(var = rownames(sens$PRCC), est = sens$PRCC$original) %>% 
  arrange(abs(est)) %>% 
  mutate(
    var = factor(var, var)
  ) %>% 
  ggplot() +
  geom_bar(aes(y = var, x = est), stat = "identity") +
  theme(axis.text.x = element_text(angle = -60, vjust = 1, hjust = 0)) +
  scale_x_continuous("PRCC to CDx sensitivity") +
  scale_y_discrete("Variable", label = lab.pars)

g_sens


g_cdx1_nb <- ggplot(post) +
  geom_point(aes(x = sens_cdx, y = spec_cdx, colour = p_scanty), alpha = 0.3) +
  scale_x_continuous("CDx sensitivity, %", labels = scales::percent) +
  scale_y_continuous("CDx specificity, %", labels = scales::percent) + 
  scale_color_continuous("Proportion TB unable to give sputum sample, %")

g_cdx1_bn <- ggplot(post) +
  geom_point(aes(x = sens_cdx_bn, y = spec_cdx_bn, colour = p_scanty), alpha = 0.3) +
  scale_x_continuous("CDx sensitivity, %", labels = scales::percent) +
  scale_y_continuous("CDx specificity, %", labels = scales::percent) + 
  scale_color_continuous("Proportion TB unable to give sputum sample, %")

g_cdx1 <- ggpubr::ggarrange(
  g_cdx1_nb + labs(title = "CDx without Bac"), 
  g_cdx1_bn + labs(title = "CDx after Bac-"), common.legend = T)

g_cdx1


g_cdx2_nb <- ggplot(post_alt) +
  geom_point(aes(x = sens_cdx, y = spec_cdx, colour = p_scanty), alpha = 0.3) +
  scale_x_continuous("CDx sensitivity, %", labels = scales::percent) +
  scale_y_continuous("CDx specificity, %", labels = scales::percent) + 
  scale_color_continuous("Proportion TB unable to give sputum sample, %")

g_cdx2_bn <- ggplot(post_alt) +
  geom_point(aes(x = sens_cdx_bn, y = spec_cdx_bn, colour = p_scanty), alpha = 0.3) +
  scale_x_continuous("CDx sensitivity, %", labels = scales::percent) +
  scale_y_continuous("CDx specificity, %", labels = scales::percent) + 
  scale_color_continuous("Proportion TB unable to give sputum sample, %")


g_cdx2 <- ggpubr::ggarrange(
  g_cdx2_nb + labs(title = "CDx without Bac"), 
  g_cdx2_bn + labs(title = "CDx after Bac-"), common.legend = T)

tab


ggsave(g_spec, filename = here::here("docs", "figs", "g_cdx_spec.png"), width = 4, height = 4)
ggsave(g_sens, filename = here::here("docs", "figs", "g_cdx_sens.png"), width = 4, height = 4)
ggsave(g_cdx, filename = here::here("docs", "figs", "g_cdx_cross.png"), width = 5, height = 5)

ggsave(g_cdx1, filename = here::here("docs", "figs", "g_cdx_cross1.png"), width = 8, height = 5)
ggsave(g_cdx2, filename = here::here("docs", "figs", "g_cdx_cross2.png"), width = 8, height = 5)



post_alt %>% 
  select(starts_with("ppv_")) %>% 
  pivot_longer(everything()) %>% 
  ggplot() +
  stat_halfeye(aes(x = value, y = name)) +
  scale_x_continuous("PPV, %", labels = scales::percent) +
  expand_limits(x = c(0, 1))



post_alt %>% 
  select(starts_with("pdx_")) %>% 
  pivot_longer(everything()) %>% 
  ggplot() +
  stat_halfeye(aes(x = value, y = name)) +
  scale_x_continuous("P(Dx), %", labels = scales::percent) +
  expand_limits(x = c(0, 1))



post_alt %>% 
  ggplot() +
  stat_halfeye(aes(x = p_scanty, fill = "2CDx"), alpha = 0.2) +
  #stat_halfeye(data=post, aes(x = p_scanty, fill = "1CDx"), alpha = 0.2) +
  scale_x_continuous("Proportion TB unable to give sputum sample, %", labels = scales::percent) +
  scale_y_continuous("Density") +
  scale_fill_discrete("")
  expand_limits(x = c(0.1, 0.6))

  
  
post_alt %>% 
  ggplot() +
  stat_halfeye(aes(x = p_scanty), alpha = 0.7) +
  geom_vline(xintercept = 0.15, linetype = 2) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  #stat_halfeye(data=post, aes(x = p_scanty, fill = "1CDx"), alpha = 0.2) +
  scale_x_continuous("Proportion TB unable to give sputum sample, %", labels = scales::percent) +
    scale_y_continuous("Density") +
  expand_limits(x = c(0.1, 0.6))

post %>% 
  select(p_scanty) %>% 
  pivot_longer(everything()) %>% 
  ggplot() +
  stat_halfeye(aes(x = value, y = name)) +
  scale_x_continuous("Fail to complete Bac or Bac-, %", labels = scales::percent) +
  expand_limits(x = c(0.1, 0.6))




post %>% 
  select(starts_with(c("sens_cdx", "spec_cdx")), p_scanty) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(
    M = median(value) * 100,
    L = quantile(value, 0.025) * 100,
    U = quantile(value, 0.975) * 100
  ) %>% 
  mutate(
    MLU = sprintf("%.0f (%.0f - %.0f)", M, L, U)
  )


post_alt %>% 
  select(starts_with(c("sens_cdx", "spec_cdx")), p_scanty) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(
    M = median(value) * 100,
    L = quantile(value, 0.025) * 100,
    U = quantile(value, 0.975) * 100
  ) %>%  
  mutate(
    MLU = sprintf("%.0f (%.0f - %.0f)", M, L, U)
  )


