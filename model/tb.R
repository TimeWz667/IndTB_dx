## ODEs -----
deriv(Y[, , ]) <- adj[k] * Y[i, j, k] + d_pop[i, j, k] + d_tb[i, j, k] + d_care[i, j, k]
dim(Y) <- c(n_tbp, n_drp, n_agp)


adj[] <- if (t > 2000) 0 else - sum(d_pop[, , i]) / sum(Y[, , i])
dim(adj) <- n_agp

deriv(CumInc[, ]) <- inc[i, j]
dim(CumInc) <- c(n_drp, n_agp)

deriv(CumMor[, ]) <- sum(death_tb[, i, j])
dim(CumMor) <- c(n_drp, n_agp)

## Initial values -----
initial(Y[, , ]) <- Y0[i, j, k]

Y0[, , ] <- user()
dim(Y0) <- c(n_tbp, n_drp, n_agp)

initial(CumInc[, ]) <- 0
initial(CumMor[, ]) <- 0

## Output -----
output(N) <- n
output(N_Age[]) <- sum(Y[, , i])
dim(N_Age) <- n_agp

output(Inc) <- sum(inc) / n * 1E5

output(Net_Pop) <- sum(d_pop)
output(Net_TB) <- sum(d_tb)
output(Net_Care) <- sum(d_care)
# output(Net_acf) <- sum(d_acf)
# output(Net_vac) <- sum(d_vac)


## Summary -----
n_tb_ut <- n_tb_a + n_tb_s + n_tb_c
n_tb_a <- sum(Y[I_Asym, , ])
n_tb_s <- sum(Y[I_Sym, , ])
n_tb_c <- sum(Y[I_ExCS, , ])
n_tb_r <- sum(Y[I_ReCS, , ])

output(PrevUt) <- n_tb_ut / n
output(PrevA) <- n_tb_ut / n
output(PrevS) <- n_tb_ut / n
output(PrevC) <- n_tb_ut / n
output(PrevR) <- n_tb_ut / n

output(PrIncRecent) <- sum(act) / sum(inc)
output(PrIncRelapse) <- (sum(rel_hi_pu) + sum(rel_hi_pr) + sum(rel_st_pu) + sum(rel_st_pr)) / sum(inc)
output(PrLat) <- (sum(Y[I_FLat, , ]) + sum(Y[I_SLat, , ]) + sum(Y[I_RHighPub, , ]) + 
                    sum(Y[I_RStPub, , ]) + sum(Y[I_RHighPri, , ]) + sum(Y[I_RStPri, , ])) / n

n <- sum(Y)

ns[, ] <- sum(Y[, i, j])
dim(ns) <- c(n_drp, n_agp)


# TB status
I_U <- 1
I_Asym <- 2
I_Sym <- 3 # Symptomatic without care-seeking
I_ExCS <- 4 # Revisit care-seeking
I_ReCS <- 5 # Returned care-seeker after LTFU from tx
I_TxPub <- 6
I_TxPriOnPub <- 7
I_TxPriOnPri <- 8
I_TxnPub <- 9 # New treatment, public
I_TxnPriOnPub <- 10  # New treatment, engaged private

I_FLat <- 11
I_SLat <- 12
I_RHighPub <- 13
I_RStPub <- 14
I_RHighPri <- 15
I_RStPri <- 16

## lengths -----
n_tbp <- 16
n_agp <- user(7)
n_drp <- 3


## population dynamics -----

tt_demo[] <- user() # data times for AIM by 0.1
dim(tt_demo) <- user() # note need this before length() use
n_tt_demo <- length(tt_demo)

r_birth[] <- user()
dim(r_birth) <- n_tt_demo
br_t <- interpolate(tt_demo, r_birth, "linear")

r_death[, ] <- user()
dim(r_death) <- c(n_tt_demo, n_agp)
dr_t[] <- interpolate(tt_demo, r_death, "linear")
dim(dr_t) <- n_agp

r_die_sym <- user(0.12)
rr_die_asym <- user(1)
r_die_asym <- rr_die_asym * r_die_sym

r_aging[, ] <- user()
dim(r_aging) <- c(n_tt_demo, n_agp)
ar_t[] <- interpolate(tt_demo, r_aging, "linear")
dim(ar_t) <- n_agp

n0[, ] <- user()
dim(n0) <- c(n_tt_demo, n_agp)
n0_t[] <- interpolate(tt_demo, n0, "linear")
dim(n0_t) <- n_agp

mr_t[] <- 50 * (n0_t[i] - sum(Y[, , i])) / n0_t[i]
dim(mr_t) <- n_agp


dr_tb[I_Asym, ] <- r_die_asym
dr_tb[I_Sym:I_ReCS, ] <- r_die_sym
dim(dr_tb) <- c(n_tbp, n_drp)

death_tb[, , ] <- dr_tb[i, j] * Y[i, j, k]
dim(death_tb) <- c(n_tbp, n_drp, n_agp)

mu_t[] <- dr_t[i] - sum(death_tb[, , i]) / sum(Y[, , i])
dim(mu_t) <- n_agp

d_pop[1, 1, 1] <- br_t * n - (mu_t[k] + dr_tb[i, j] + ar_t[k] - mr_t[k]) * Y[i, j, k]
d_pop[2:n_tbp, 1:3, 1] <-  - (mu_t[k] + dr_tb[i, j] + ar_t[k] - mr_t[k]) * Y[i, j, k]
d_pop[1, 2:3, 1] <-        - (mu_t[k] + dr_tb[i, j] + ar_t[k] - mr_t[k]) * Y[i, j, k]

d_pop[, , 2:n_agp] <- ar_t[k - 1] * Y[i, j, k - 1] - (mu_t[k] + dr_tb[i, j] + ar_t[k] - mr_t[k]) * Y[i, j, k]
dim(d_pop) <- c(n_tbp, n_drp, n_agp)



## TB -----

beta <- user(10)
rr_beta_dr <- user(0.6)
rr_beta_nr <- rr_beta_dr

rr_sus_slat <- user(0.2)

r_lat <- user(0.5)
r_stab <- user(0.667)

p_primary <- user(0.1)
r_act <- r_lat * p_primary / (1 - p_primary)


r_clear <- user(0.03)

r_react <- user(0.001) # Activation of Latent TB
r_relapse <- user(0.001) # relapse after stabilisation
r_relapse_te <- 0.1487 * 1.5 / 2 / 0.6

r_relapse_new_tx <- user(0.1)
pr_txe_new_pu <- sum(txe[4, , ]) / (sum(txe[1, , ]) + sum(txe[4, , ]))
pr_txe_new_pr <- sum(txe[5, , ]) / (sum(txe[2, , ]) + sum(txe[3, , ]) + sum(txe[5, , ]))

r_rel_hi_pu <- min(pr_txe_new_pu * r_relapse_new_tx + (1 - pr_txe_new_pu) * r_relapse_te, r_relapse_te)
r_rel_hi_pr <- min(pr_txe_new_pr * r_relapse_new_tx + (1 - pr_txe_new_pr) * r_relapse_te, r_relapse_te)


r_sc <- user(0.2)
r_onset <- user(1)

## Transitions and natural history
#### Intervention vaccination
tt_intv_vac[] <- user() # uptakes for mass screening intervention
dim(tt_intv_vac) <- user()

uptake_vac[] <- user()
dim(uptake_vac) <-user()
uptake_vac_t <- interpolate(tt_intv_vac, uptake_vac, "linear")

cover_vac <- user(1)

prot_inf <- user(0)
prot_prog <- user(0)
prot_rel <- user(0)

k_inf <- (1 - uptake_vac_t) + uptake_vac_t * (1 - prot_inf * cover_vac)




### TB infection
irr_25 <- user(1)
irr_35 <- user(1)
irr_45 <- user(1)
irr_55 <- user(1)
irr_65 <- user(1)


foi[1] <- beta              * (sum(Y[I_Asym, 1, ]) + sum(Y[I_Sym, 1, ]) + sum(Y[I_ExCS, 1, ]) + sum(Y[I_ReCS, 1, ])) / n
foi[2] <- beta * rr_beta_dr * (sum(Y[I_Asym, 2, ]) + sum(Y[I_Sym, 2, ]) + sum(Y[I_ExCS, 2, ]) + sum(Y[I_ReCS, 2, ])) / n
foi[3] <- beta * rr_beta_nr * (sum(Y[I_Asym, 3, ]) + sum(Y[I_Sym, 3, ]) + sum(Y[I_ExCS, 3, ]) + sum(Y[I_ReCS, 3, ])) / n
dim(foi) <- 3


inf_u[, 2, ] <-  k_inf * foi[k] * Y[I_U, i, j]
inf_u[, 3, ] <- irr_25 * foi[k] * Y[I_U, i, j]
inf_u[, 4, ] <- irr_35 * foi[k] * Y[I_U, i, j]
inf_u[, 5, ] <- irr_45 * foi[k] * Y[I_U, i, j]
inf_u[, 6, ] <- irr_55 * foi[k] * Y[I_U, i, j]
inf_u[, 7, ] <- irr_65 * foi[k] * Y[I_U, i, j]
dim(inf_u) <- c(n_drp, n_agp, n_drp)

inf_sl[, 2, ] <-  k_inf * foi[k] * rr_sus_slat * Y[I_SLat, i, j]
inf_sl[, 3, ] <- irr_25 * foi[k] * rr_sus_slat * Y[I_SLat, i, j]
inf_sl[, 4, ] <- irr_35 * foi[k] * rr_sus_slat * Y[I_SLat, i, j]
inf_sl[, 5, ] <- irr_45 * foi[k] * rr_sus_slat * Y[I_SLat, i, j]
inf_sl[, 6, ] <- irr_55 * foi[k] * rr_sus_slat * Y[I_SLat, i, j]
inf_sl[, 7, ] <- irr_65 * foi[k] * rr_sus_slat * Y[I_SLat, i, j]
dim(inf_sl) <- c(n_drp, n_agp, n_drp)

inf_ru[, 2, ] <-  k_inf * foi[k] * rr_sus_slat * Y[I_RStPub, i, j]
inf_ru[, 3, ] <- irr_25 * foi[k] * rr_sus_slat * Y[I_RStPub, i, j]
inf_ru[, 4, ] <- irr_35 * foi[k] * rr_sus_slat * Y[I_RStPub, i, j]
inf_ru[, 5, ] <- irr_45 * foi[k] * rr_sus_slat * Y[I_RStPub, i, j]
inf_ru[, 6, ] <- irr_55 * foi[k] * rr_sus_slat * Y[I_RStPub, i, j]
inf_ru[, 7, ] <- irr_65 * foi[k] * rr_sus_slat * Y[I_RStPub, i, j]
dim(inf_ru) <- c(n_drp, n_agp, n_drp)

inf_ri[, 2, ] <-  k_inf * foi[k] * rr_sus_slat * Y[I_RStPri, i, j]
inf_ri[, 3, ] <- irr_25 * foi[k] * rr_sus_slat * Y[I_RStPri, i, j]
inf_ri[, 4, ] <- irr_35 * foi[k] * rr_sus_slat * Y[I_RStPri, i, j]
inf_ri[, 5, ] <- irr_45 * foi[k] * rr_sus_slat * Y[I_RStPri, i, j]
inf_ri[, 6, ] <- irr_55 * foi[k] * rr_sus_slat * Y[I_RStPri, i, j]
inf_ri[, 7, ] <- irr_65 * foi[k] * rr_sus_slat * Y[I_RStPri, i, j]
dim(inf_ri) <- c(n_drp, n_agp, n_drp)


### TB progression
act[, ] <- r_act * Y[I_FLat, i, j]
dim(act) <- c(n_drp, n_agp)

react[, ] <- r_react * Y[I_SLat, i, j]
dim(react) <- c(n_drp, n_agp)

rel_hi_pu[, ] <- r_rel_hi_pu * Y[I_RHighPub, i, j]
dim(rel_hi_pu) <- c(n_drp, n_agp)

rel_hi_pr[, ] <- r_rel_hi_pr * Y[I_RHighPri, i, j]
dim(rel_hi_pr) <- c(n_drp, n_agp)

rel_st_pu[, ] <- r_relapse * Y[I_RStPub, i, j]
dim(rel_st_pu) <- c(n_drp, n_agp)

rel_st_pr[, ] <- r_relapse * Y[I_RStPri, i, j]
dim(rel_st_pr) <- c(n_drp, n_agp)

# inc[, ] <- incr * ns[i, j]
inc[, ] <- act[i, j] + react[i, j] + rel_hi_pu[i, j] + rel_hi_pr[i, j] + rel_st_pu[i, j] + rel_st_pr[i, j]
dim(inc) <- c(n_drp, n_agp)


d_tb[I_U, , ] <- - sum(inf_u[j, k, ]) + r_clear * (Y[I_RStPub, j, k] + Y[I_RStPri, j, k] + Y[I_SLat, j, k])
d_tb[I_Asym, , ] <- inc[j, k] - (r_onset + r_sc) * Y[I_Asym, j, k]
d_tb[I_Sym, , ] <- r_onset * Y[I_Asym, j, k] - r_sc * Y[I_Sym, j, k]
d_tb[I_ExCS, , ] <- - r_sc * Y[I_ExCS, j, k]
d_tb[I_ReCS, , ] <- - r_sc * Y[I_ReCS, j, k]

d_tb[I_RHighPub, , ] <- - rel_hi_pu[j, k] - r_stab * Y[i, j, k]
d_tb[I_RHighPri, , ] <- - rel_hi_pr[j, k] - r_stab * Y[i, j, k]
d_tb[I_RStPub, , ] <- r_stab * Y[I_RHighPub, j, k] - rel_st_pu[j, k] - r_clear * Y[I_RStPub, j, k] - sum(inf_ru[j, k, ])
d_tb[I_RStPri, , ] <- r_stab * Y[I_RHighPri, j, k] - rel_st_pr[j, k] - r_clear * Y[I_RStPri, j, k] - sum(inf_ri[j, k, ])

d_tb[I_FLat, , ] <- sum(inf_u[, k, j]) + sum(inf_sl[, k, j]) + sum(inf_ru[, k, j]) + sum(inf_ri[, k, j]) - r_lat * Y[I_FLat, j, k] - act[j, k]
d_tb[I_SLat, , ] <- r_lat * Y[I_FLat, j, k] - react[j, k] - r_clear * Y[I_SLat, j, k] - sum(inf_sl[j, k, ]) + 
                    r_sc * (Y[I_Asym, j, k] + Y[I_Sym, j, k] + Y[I_ExCS, j, k] + Y[I_ReCS, j, k])

dim(d_tb) <- c(n_tbp, n_drp, n_agp)



### Intervention Diagnosis
r_csi0 <- user(3)
r_recsi0 <- user(6)

t0_decline <- 2005
rt_cs <- user(0.05)

k_cs <- if (t < t0_decline) 1 else exp(rt_cs * (t - t0_decline))
r_csi <- r_csi0 * k_cs
r_recsi <- r_recsi0 * k_cs


tt_intv_dx[] <- user() # uptakes for treatment intervention
dim(tt_intv_dx) <- user()

uptake_dx[] <- user()
dim(uptake_dx) <- user()
uptake_dx_t <- interpolate(tt_intv_dx, uptake_dx, "linear")


p_ent[] <- user()
dim(p_ent) <- 3

p_itt <- user(0.5)

p_dx[] <- user()
dim(p_dx) <- 3

p_txi[] <- user()
dim(p_txi) <- 3

p_pri_on_pub <- user()

p_cure[, ] <- user()
dim(p_cure) <- c(3, 3)


p_det[] <- p_ent[i] * p_itt * p_dx[i]
dim(p_det) <- 3 # sectors

p_fn <- 1 - sum(p_det)

p_txs_soc[1, ] <- p_det[1] * (1 - uptake_tx_t) * p_txi[1] *                      p_cure[1, j]
p_txs_soc[2, ] <- p_det[2] * (1 - uptake_tx_t) * p_txi[2] *      p_pri_on_pub *  p_cure[2, j] + 
                  p_det[3] *                     p_txi[3] * (1 - p_pri_on_pub) * p_cure[3, j]
p_txs_soc[3, ] <- p_det[3] * p_txi[3] * p_cure[3, j]
dim(p_txs_soc) <- c(3, n_drp)

tx_dur[] <- user()
dim(tx_dur) <- 3

### Intervention Treatment
tt_intv_tx[] <- user() # uptakes for treatment intervention
dim(tt_intv_tx) <- user()

uptake_tx[] <- user()
dim(uptake_tx) <- user()
uptake_tx_t <- interpolate(tt_intv_tx, uptake_tx, "linear")

p_txi_new[] <- user()
dim(p_txi_new) <- 3

p_cure_new[, ] <- user()
dim(p_cure_new) <- c(3, 3)

p_txs_new[1, ] <- p_det[1] * uptake_tx_t * p_txi_new[1] *                p_cure_new[1, j]
p_txs_new[2, ] <- p_det[2] * uptake_tx_t * p_txi_new[2] * p_pri_on_pub * p_cure_new[2, j]
p_txs_new[3, ] <- 0 # no new treatment in non-engaged private

dim(p_txs_new) <- c(3, n_drp)

tx_new_dur <- user(1 / 12)

p_txs[, ] <- p_txs_soc[i, j] + p_txs_new[i, j]
dim(p_txs) <- c(3, n_drp)


### Intervention ACF/Mass-screening
tt_intv_ms[] <- user() # uptakes for mass screening intervention
dim(tt_intv_ms) <- user()

uptake_ms[] <- user()
dim(uptake_ms) <- user()
uptake_ms_t <- interpolate(tt_intv_ms, uptake_ms, "linear")

cover_acf <- user(0.3)
sens_acf <- user(0.7)

acf_with_new_tx <- user(0)

r_acf <- cover_acf * sens_acf * uptake_ms_t

r_acf_soc <- r_acf * (1 - acf_with_new_tx) # ACF with soc treatment
r_acf_new <- r_acf * acf_with_new_tx # ACF with intervening treatment


## TB Care
fn0[, ] <- r_csi * p_fn * Y[I_Sym, i, j]
dim(fn0) <- c(n_drp, n_agp)

ts_soc0[, , ] <- r_csi * p_txs_soc[i, j] * Y[I_Sym, j, k]
dim(ts_soc0) <- c(3, n_drp, n_agp)

ts_new0[, , ] <- r_csi * p_txs_new[i, j] * Y[I_Sym, j, k]
dim(ts_new0) <- c(3, n_drp, n_agp)

tl0[, ] <- r_csi * Y[I_Sym, i, j] - fn0[i, j] - sum(ts_soc0[, i, j]) - sum(ts_new0[, i, j])
dim(tl0) <- c(n_drp, n_agp)

fn1[, ] <- r_recsi * p_fn * Y[I_ExCS, i, j]
dim(fn1) <- c(n_drp, n_agp)

ts_soc1[, , ] <- r_recsi * p_txs_soc[i, j] * Y[I_ExCS, j, k]
dim(ts_soc1) <- c(3, n_drp, n_agp)

ts_new1[, , ] <- r_recsi * p_txs_new[i, j] * Y[I_ExCS, j, k]
dim(ts_new1) <- c(3, n_drp, n_agp)

tl1[, ] <- r_recsi * Y[I_ExCS, i, j] - fn1[i, j] - sum(ts_soc1[, i, j]) - sum(ts_new1[, i, j])
dim(tl1) <- c(n_drp, n_agp)

fn2[, ] <- r_recsi * p_fn * Y[I_ReCS, i, j]
dim(fn2) <- c(n_drp, n_agp)

ts_soc2[, , ] <- r_recsi * p_txs_soc[i, j] * Y[I_ReCS, j, k]
dim(ts_soc2) <- c(3, n_drp, n_agp)

ts_new2[, , ] <- r_recsi * p_txs_new[i, j] * Y[I_ReCS, j, k]
dim(ts_new2) <- c(3, n_drp, n_agp)

tl2[, ] <- r_recsi * Y[I_ReCS, i, j] - fn2[i, j] - sum(ts_soc2[, i, j]) - sum(ts_new2[, i, j])
dim(tl2) <- c(n_drp, n_agp)

txe[1, , ] <- (1 / tx_dur[1]) * Y[I_TxPub, j, k]
txe[2, , ] <- (1 / tx_dur[2]) * Y[I_TxPriOnPub, j, k]
txe[3, , ] <- (1 / tx_dur[3]) * Y[I_TxPriOnPri, j, k]
txe[4, , ] <- (1 / tx_new_dur) * Y[I_TxnPub, j, k]
txe[5, , ] <- (1 / tx_new_dur) * Y[I_TxnPriOnPub, j, k]

dim(txe) <- c(5, n_drp, n_agp)

d_care[I_Asym, , ] <- - r_acf * Y[I_Asym, j, k]
d_care[I_Sym, , ] <- - (r_csi + r_acf) * Y[I_Sym, j, k]
d_care[I_ExCS, , ] <- fn0[j, k] + fn1[j, k] - (r_recsi + r_acf) * Y[I_ExCS, j, k]
d_care[I_ReCS, , ] <- fn2[j, k] + tl0[j, k] + tl1[j, k] + tl2[j, k] - (r_recsi + r_acf) * Y[I_ReCS, j, k]

d_care[I_TxPub, , ]      <- ts_soc0[1, j, k] + ts_soc1[1, j, k] + ts_soc2[1, j, k] - txe[1, j, k] + 
                             r_acf_soc * (Y[I_Asym, j, k] + Y[I_Sym, j, k] + Y[I_ExCS, j, k] + Y[I_ReCS, j, k])
d_care[I_TxPriOnPub, , ] <- ts_soc0[2, j, k] + ts_soc1[2, j, k] + ts_soc2[2, j, k] - txe[2, j, k]
d_care[I_TxPriOnPri, , ] <- ts_soc0[3, j, k] + ts_soc1[3, j, k] + ts_soc2[3, j, k] - txe[3, j, k]

d_care[I_TxnPub, , ]      <- ts_new0[1, j, k] + ts_new1[1, j, k] + ts_new2[1, j, k] - txe[4, j, k] + 
                             r_acf_new * (Y[I_Asym, j, k] + Y[I_Sym, j, k] + Y[I_ExCS, j, k] + Y[I_ReCS, j, k])
d_care[I_TxnPriOnPub, , ] <- ts_new0[2, j, k] + ts_new1[2, j, k] + ts_new2[2, j, k] - txe[5, j, k]

d_care[I_RHighPub, , ] <- txe[1, j, k] + txe[4, j, k]
d_care[I_RHighPri, , ] <- txe[2, j, k] + txe[3, j, k] + txe[5, j, k]

dim(d_care) <- c(n_tbp, n_drp, n_agp)
