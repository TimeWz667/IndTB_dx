library(tidyverse)



# ## From review
# delay_sys <- c(31.0 + 2.5, 24.5 + 1.9, 35.4 + 3.6)
n_vis <- c(2.7, 1.9, 12.3)

opt_vis <- nlminb(c(1, 1), function(x) {
  qs <- qlnorm(c(0.5, 0.25, 0.75), x[1], x[2])
  sum((qs / n_vis - 1) ** 2)
}, lower = 0)
opt_vis





## Load data
for (folder in c("bac_cdx_sector_2021", "bac_cdx_sector_2022")) {
  js <- jsonlite::read_json(here::here("pars", "pars_" + glue::as_glue(folder) + ".json"))
  prev <- js$prev
  
  post <- read_csv(here::here("out", folder, "post.csv"))
  
  red <- post %>% 
    mutate(
      r_tb = r_test * p_tb,
      p_ent_pub = p_csi_pub,
      p_ent_eng = (1 - p_csi_pub) * p_csi_ppm,
      p_ent_pri = (1 - p_csi_pub) * (1 - p_csi_ppm),
      p_itt_pub = 1,
      p_itt_eng = rp_ent,
      p_itt_pri = rp_ent,
      
      p_dx_pub = p_tp_bac_pub + p_tp_cdx_pub,
      p_dx_eng = p_tp_bac_eng + p_tp_cdx_eng,
      p_dx_pri = p_tp_bac_pub + p_tp_cdx_pub,
      r_txi_pub = r_tb_pub * p_dx_pub * p_txi_pub,
      r_txi_eng = r_tb_eng * p_dx_eng * p_txi_eng,
      r_txi_pri = r_tb_pri * p_dx_pri * p_txi_pri,
      r_txi = r_txi_pub + r_txi_eng + r_txi_pri,
      cas_itt = p_ent_pub * p_itt_pub + p_ent_eng * p_itt_eng + p_ent_pri * p_itt_pri,
      cas_dx = p_ent_pub * p_itt_pub * p_dx_pub + p_ent_eng * p_itt_eng * p_dx_eng + p_ent_pri * p_itt_pri * p_dx_pri,
      cas_txi = p_ent_pub * p_itt_pub * p_dx_pub * p_txi_pub + 
        p_ent_eng * p_itt_eng * p_dx_eng * p_txi_eng + 
        p_ent_pri * p_itt_pri * p_dx_pri * p_txi_pri,
      cas_txi_full = r_txi / r_tb,
      prev_a = prev$PrevAsym,
      prev_s = prev$PrevSym,
      prev_c = prev$PrevExCS,
      r_die_a = 0.127,
      r_die_s = 0.127,
      r_sc = 0.2,
      r_det = r_txi / prev_c,
      r_csi = (r_det + r_sc + r_die_s) * prev_c / prev_s,
      r_onset = (r_csi + r_sc + r_die_s) * prev_s / prev_a,
      inc = (r_onset + r_sc + r_die_a) * prev_a
    )
  
  cas_txi <- red %>% summarise(N = mean(cas_txi_full)) %>% unlist() %>% c
  k_itt <- 1 / unname(cas_txi * n_vis[1])
  
  
  reformed <- red %>% 
    mutate(
      max_rat0 = r_txi / (prev_s * r_csi * cas_txi),
      rat0 = pmin(max_rat0, 1),
      # k_itt = rbinom(n(), prob = 0.21, size = 250) / 250 / rp_ent,
      # k_itt = rbinom(n(), prob = 0.74, size = 88) / 88,
      p_itt0_pub = p_itt_pub * rat0 * k_itt,
      p_itt0_eng = p_itt_eng * rat0 * k_itt,
      p_itt0_pri = p_itt_pri * rat0 * k_itt,
      p_itt1_pub = p_itt_pub * k_itt,
      p_itt1_eng = p_itt_eng * k_itt,
      p_itt1_pri = p_itt_pri * k_itt,
      p0 = cas_txi * k_itt * rat0,
      p1 = cas_txi * k_itt,
      det0 = prev_s * r_csi * p0,
      det1 = r_txi - det0,
      r_recsi = det1 / (prev_c * p1)
    ) %>% 
    select(-c(det0, det1), -starts_with(c("p_itt_", "prev_")))
  
  
  tab <- reformed %>% 
    pivot_longer(c(starts_with("r_"), starts_with("p_"), starts_with("alg"), p0, p1)) %>% 
    group_by(name) %>% 
    summarise(
      mu = mean(value),
      std = sd(value),
      m = median(value),
      l = quantile(value, 0.25),
      u = quantile(value, 0.75)
    )

  tab
  
  write_csv(reformed, file = here::here("out", folder, "post_filled.csv"))
  write_csv(tab, file = here::here("out", folder, "summary_filled.csv"))
  
  js_re <- list(
    pars = reformed,
    prev = js$prev,
    txo = js$txo
  )
  
  jsonlite::write_json(js_re, here::here("pars", "pars_" + glue::as_glue(folder) + "_re.json"), digits = 8, auto_unbox = T)
  
}

