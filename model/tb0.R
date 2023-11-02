## ODEs -----
deriv(Y[, , , ]) <- d_pop[i, j, k, l]+ d_hiv[i, j, k, l] + d_tb[i, j, k, l]
dim(Y) <- c(n_agp, 2, 3, n_tbp)


## Initial values -----
initial(Y[, , , ]) <- Y0[i, j, k, l]

Y0[, , , ] <- user()
dim(Y0) <- c(n_agp, 2, 3, n_tbp)


# 12 states for tuberculosis:
# 1: sus, 2: fast latent, 3: slow latent, 
# 4: sub-clinical tb, 5: sm-tb pre-care-seeking, 6: sm+tb pre-care-seeking
# 7: tb on treating
# 8: fast latent with treatment history, 9: recovered
# 10: re sub-clinical tb, 11: re sm-tb pre-care-seeking, 12: re sm+tb pre-care-seeking


## Output -----
output(N_HIV) <- poph
output(N_ART) <- popart
output(N_1549) <- pop1549
output(N_All) <- n_all
output(N_F) <- n_f
output(N_M) <- n_m
output(N_TB) <- n_tb
output(N_TB_Sub) <- n_tb_sub
output(N_TB_Sn) <- n_tb_sn
output(N_TB_Sp) <- n_tb_sp
output(N_TB_HIV) <- n_tb_hiv

output(Inc) <- sum(inc) / n_all * 1E5
output(Noti) <- (sum(noti_sp) + sum(noti_sn)) / n_all * 1E5

output(Prv_F) <- n_tb_f * 1E5
output(Inc_F) <- sum(inc[, 1, ]) / n_f * 1E5
output(Noti_F) <- (sum(noti_sp[, 1, ]) + sum(noti_sn[, 1, ])) / n_f * 1E5
output(Noti_F15U) <- (sum(noti_sp[4:n_agp, 1, ]) + sum(noti_sn[4:n_agp, 1, ])) / n_f15u * 1E5
output(Noti_F15U_Sp) <- sum(noti_sp[4:n_agp, 1, ]) / n_f15u * 1E5


output(Prv_M) <- n_tb_m * 1E5
output(Inc_M) <- sum(inc[, 2, ]) / n_m * 1E5
output(Noti_M) <- (sum(noti_sp[, 2, ]) + sum(noti_sn[, 2, ])) / n_m * 1E5
output(Noti_M15U) <- (sum(noti_sp[4:n_agp, 2, ]) + sum(noti_sn[4:n_agp, 2, ])) / n_m15u * 1E5
output(Noti_M15U_Sp) <- sum(noti_sp[4:n_agp, 2, ]) / n_m15u * 1E5

output(Det_Sn_F) <- r_det_sn[1]
output(Det_Sn_M) <- r_det_sn[2]
output(Det_Sp_F) <- r_det_sp[1]
output(Det_Sp_M) <- r_det_sp[2]

output(Net_TB) <- sum(d_tb)
output(Net_HIV) <- sum(d_hiv)
output(Net_Pop) <- sum(d_pop)
output(Ratio_Recent) <- pr_recent
output(Pr_Lat) <- (sum(Y[, , , 2:3]) + sum(Y[, , , 8:9])) / n_all

## Summary -----
pop1549 <- sum(Y[4:10, , , ])
poph1549 <- sum(Y[4:10, , 2:3, ])
poph <- sum(Y[4:n_agp, , 2:3, ])
popart <- sum(Y[4:n_agp, , 3, ])

n_tb <- n_tb_sub + n_tb_sn + n_tb_sp
n_tb_sub <- sum(Y[, , , 4]) + sum(Y[, , , 10])
n_tb_sn <- sum(Y[, , , 5]) + sum(Y[, , , 11])
n_tb_sp <- sum(Y[, , , 6]) + sum(Y[, , , 12])
n_tb_f <- sum(Y[, 1, , 4:6]) + sum(Y[, 1, , 10:12])
n_tb_m <- sum(Y[, 2, , 4:6]) + sum(Y[, 2, , 10:12])
n_tb_hiv <- sum(Y[, , 2:3, 4:6]) + sum(Y[, , 2:3, 10:12])
n_f <- sum(Y[, 1, , ])
n_m <- sum(Y[, 2, , ])
n_f15u <- sum(Y[4:n_agp, 1, , ])
n_m15u <- sum(Y[4:n_agp, 2, , ])
n_all <- sum(Y)

pr_recent <- r_act * (sum(Y[, , , 2]) + sum(Y[, , , 8])) / (sum(inc[, , ]) + sum(re_inc[, , ]))


## dims -----
ageing[] <- user()
dim(ageing) <- n_agp

ttp[] <- user() # data times for AIM by 0.1
dim(ttp) <- user() # note need this before length() use
tt1[] <- user() # data times for AIM by 1
dim(tt1) <- user()

tscale <- user()

## lengths -----
n_agp <- user()
n_tbp <- 12
n_tt1 <- length(tt1)
n_ttp <- length(ttp)

## birth -----
br[, ] <- user()
dim(br) <- c(n_tt1, 2)

br_t[] <- interpolate(tt1, br, "linear")
dim(br_t) <- 2

## death + migration -----
omega[, , ] <- user()
dim(omega) <- c(n_tt1, n_agp, 2)

mu_h_f[, ] <- user()
mu_h_m[, ] <- user()
dim(mu_h_f) <- c(n_ttp, n_agp)
dim(mu_h_m) <- c(n_ttp, n_agp)
mu_h[, , 1] <- mu_h_f[i, j]
mu_h[, , 2] <- mu_h_m[i, j]
dim(mu_h) <- c(n_ttp, n_agp, 2)

mu_a_f[, ] <- user()
mu_a_m[, ] <- user()
dim(mu_a_f) <- c(n_ttp, n_agp)
dim(mu_a_m) <- c(n_ttp, n_agp)
mu_a[, , 1] <- mu_a_f[i, j]
mu_a[, , 2] <- mu_a_m[i, j]
dim(mu_a) <- c(n_ttp, n_agp, 2)

omega_t[, ] <- interpolate(tt1, omega, "linear")
dim(omega_t) <- c(n_agp, 2)

mu_h_t[, ] <- interpolate(ttp, mu_h, "linear")
dim(mu_h_t) <- c(n_agp, 2)

mu_a_t[, ] <- interpolate(ttp, mu_a, "linear")
dim(mu_a_t) <- c(n_agp, 2)

mu_tb[1:4] <- 0
mu_tb[5] <- r_tbmu_sn_untr
mu_tb[6] <- r_tbmu_sp_untr
mu_tb[7] <- r_tbmu_ontr
mu_tb[8:10] <- 0
mu_tb[11] <- r_tbmu_sn_untr
mu_tb[12] <- r_tbmu_sp_untr

dim(mu_tb) <- n_tbp

mu_hiv[, , 1] <- 0
mu_hiv[, , 2] <- mu_h_t[i, j]
mu_hiv[, , 3] <- mu_a_t[i, j]
dim(mu_hiv) <- c(n_agp, 2, 3) 

mu_x_t[, ] <- (omega_t[i, j] * sum(Y[i, j, , ]) 
               - mu_h_t[i, j] * sum(Y[i, j, 2, ]) 
               - mu_a_t[i, j] * sum(Y[i, j, 3, ])
               - r_tbmu_sn_untr * (sum(Y[i, j, , 5]) + sum(Y[i, j, , 11]))
               - r_tbmu_sp_untr * (sum(Y[i, j, , 6]) + sum(Y[i, j, , 12]))
               - r_tbmu_ontr * sum(Y[i, j, , 7])
) / (sum(Y[i, j, , ]))
dim(mu_x_t) <- c(n_agp, 2)  

mu_t[, , , ] <- mu_x_t[i, j] + mu_hiv[i, j, k] + mu_tb[l]
dim(mu_t) <- c(n_agp, 2, 3, n_tbp)


d_pop[1, , 1, 1] <- br_t[j] * sum(Y[, j, , ]) - (mu_t[1, j, k, l] + ageing[i]) * Y[1, j, k, l]
d_pop[1, , 1, 2:12] <- 0 - (mu_t[1, j, k, l] + ageing[i]) * Y[1, j, k, l]
d_pop[1, , 2:3, ] <- 0 - (mu_t[1, j, k, l] + ageing[i]) * Y[1, j, k, l]
d_pop[2:n_agp, , , ] <- ageing[i - 1] * Y[i - 1, j, k, l] - (mu_t[i, j, k, l] + ageing[i]) * Y[i, j, k, l]
dim(d_pop) <- c(n_agp, 2, 3, n_tbp)


## HIV -----

hit[] <- user()                        #HIV/ART targets
hat[] <- user()
hsr[] <- user()
dim(hit) <- n_tt1
dim(hat) <- n_tt1
dim(hsr) <- n_tt1

hit_t <- interpolate(tt1, hit, "linear")
hsr_t <- interpolate(tt1, hsr, "linear")


hiv_age[, ] <- user()                 #HIV by age/sex data
dim(hiv_age) <- c(n_agp, 2)
n_hiv_age[, ] <- hiv_age[i, j] * sum(Y[i, j, , ])
dim(n_hiv_age) <- c(n_agp, 2)

hivi <- 1 * max((hit_t * pop1549 - poph1549), 0) / tscale #pursuit
hfz[, 1] <- hivi * (hsr_t) * n_hiv_age[i, 1] / (sum(n_hiv_age[, 1]) + 1E-15) / (sum(Y[i, j, 1, ]) + 1E-15)
hfz[, 2] <- hivi * (1 - hsr_t) * n_hiv_age[i, 2] / (sum(n_hiv_age[, 2]) + 1E-15) / (sum(Y[i, j, 1, ]) + 1E-15)
dim(hfz) <- c(n_agp, 2)


## ART initialisation
hat_t <- interpolate(tt1, hat, "linear")
arti <- max(hat_t * poph - popart, 0) / tscale #pursuit

rra[, ] <- user()
dim(rra) <- c(n_agp, 2)

rraN[, ] <- rra[i, j] * sum(Y[i, j, 2, ])
dim(rraN) <- c(n_agp, 2)

raz[, ] <- arti * rraN[i, j] / (sum(rraN) + 1E-15) / (sum(Y[i, j, 2, ]) + 1E-15)
dim(raz) <- c(n_agp, 2)

d_hiv[, , 1, ] <- - hfz[i, j] * Y[i, j, 1, l] / sum(Y[i, j, 1, ])
d_hiv[, , 2, ] <- hfz[i, j] * Y[i, j, 1, l] / sum(Y[i, j, 1, ]) - raz[i, j] * Y[i, j, 2, l]
d_hiv[, , 3, ] <- raz[i, j] * Y[i, j, 2, l]

dim(d_hiv) <- c(n_agp, 2, 3, n_tbp)


## TB -----

### Parameters -----
beta <- user(20)
rr_beta_t <- user(0.01)
beta_lower <- user(1)
beta_sp <- max(beta * exp( - rr_beta_t * (t - 1970)), beta_lower)
beta_sn <- beta_sp * 0.22
rr_male <- user(2)

p_im <- user(0.79)

r_lat <- user(0.2)

r_act <- user(0.07)
r_ract <- r_act / 50 # user(0.0007)
r_rel <- r_act / 50 # user(0.0003)

rr_det <- user(0.1) # from WHO data (case detection rate) + intervention
y_mid <- user(2000)

r_seek_sn_f <- user(0.2)
r_seek_sn_m <- user(0.2)
r_seek_sp_f <- user(0.5319149)
r_seek_sp_m <- user(0.3610108)


r_det_sn[1] <- r_seek_sn_f / (1 + exp(- rr_det * (t - y_mid)))
r_det_sn[2] <- r_seek_sn_m / (1 + exp(- rr_det * (t - y_mid)))
dim(r_det_sn) <- 2
r_det_sp[1] <- r_seek_sp_f / (1 + exp(- rr_det * (t - y_mid)))
r_det_sp[2] <- r_seek_sp_m / (1 + exp(- rr_det * (t - y_mid)))
dim(r_det_sp) <- 2

r_sym <- user(2)
p_sp <- user(0.3)
p_sn <- 1 - p_sp
r_trans <- user(0.373) # qexp(pexp(0.006) * 52)

r_rec <- user(0.5)
r_cure <- user(0.15) # from WHO data (treatment successful) + intervention
r_loss <- user(0.2)
r_tbmu_sp_untr <- user(0.127) # qexp(pexp(0.01) * 12)
r_tbmu_sn_untr <- user(0.024) # qexp(pexp(0.002) * 12)
r_tbmu_ontr <- user(0.0157) # qexp(pexp(0.0013) * 12)


### Processes -----

foi0 <- (beta_sp * n_tb_sp + beta_sn * (n_tb_sn + n_tb_sub)) / sum(Y)
foi[1] <-  foi0
foi[2] <-  rr_male * foi0
re_foi[] <- (1 - p_im) * foi[i]
dim(foi) <- 2
dim(re_foi) <- 2

inc[, , ] <- r_act * Y[i, j, k, 2] + r_ract * Y[i, j, k, 3] 
re_inc[, , ] <- r_act * Y[i, j, k, 8] + r_rel * Y[i, j, k, 9]
noti_sn[, , ] <- r_det_sn[j] * (Y[i, j, k, 5] + Y[i, j, k, 11])
noti_sp[, , ] <- r_det_sp[j] * (Y[i, j, k, 6] + Y[i, j, k, 12])
# found_sub[, , ] <- r_det * (Y[i, j, k, 4] + Y[i, j, k, 10])

dim(inc) <- c(n_agp, 2, 3)
dim(re_inc) <- c(n_agp, 2, 3)
dim(noti_sn) <- c(n_agp, 2, 3)
dim(noti_sp) <- c(n_agp, 2, 3)


d_tb[, , , 1] <- - foi[j] * Y[i, j, k, 1]

d_tb[, , , 2] <- foi[j] * Y[i, j, k, 1] + re_foi[j] * Y[i, j, k, 3] - (r_lat + r_act) * Y[i, j, k, 2]
d_tb[, , , 3] <- r_lat * Y[i, j, k, 2] + r_rec * sum(Y[i, j, k, 4:6]) - (re_foi[j] + r_ract) * Y[i, j, k, 3]

d_tb[, , , 4] <- inc[i, j, k] - r_sym * Y[i, j, k, 4] - r_rec * Y[i, j, k, 4]
d_tb[, , , 5] <- r_sym * p_sn * Y[i, j, k, 4] - (r_det_sn[j] + r_trans + r_rec) * Y[i, j, k, 5]
d_tb[, , , 6] <- r_sym * p_sp * Y[i, j, k, 4] + r_trans * Y[i, j, k, 5] - (r_det_sp[j] + r_rec) * Y[i, j, k, 6]

d_tb[, , , 7] <- r_det_sn[j] * (Y[i, j, k, 5] + Y[i, j, k, 11]) + 
  r_det_sp[j] * (Y[i, j, k, 6] + Y[i, j, k, 12]) - (r_cure + r_rec + r_loss) * Y[i, j, k, 7]

d_tb[, , , 8] <- re_foi[j] * Y[i, j, k, 9] - (r_lat + r_act) * Y[i, j, k, 8] 
d_tb[, , , 9] <- r_lat * Y[i, j, k, 8] +  r_rec * sum(Y[i, j, k, 10:12]) + (r_cure + r_rec) * Y[i, j, k, 7] - (re_foi[j] + r_rel) * Y[i, j, k, 9]

d_tb[, , , 10] <- re_inc[i, j, k] + r_loss * Y[i, j, k, 7] - r_sym * Y[i, j, k, 10] - r_rec * Y[i, j, k, 10]
d_tb[, , , 11] <- r_sym * p_sn * Y[i, j, k, 10] - (r_det_sn[j] + r_trans + r_rec) * Y[i, j, k, 11]
d_tb[, , , 12] <- r_sym * p_sp * Y[i, j, k, 10] + r_trans * Y[i, j, k, 11] - (r_det_sp[j] + r_rec) * Y[i, j, k, 12]


dim(d_tb) <- c(n_agp, 2, 3, n_tbp)