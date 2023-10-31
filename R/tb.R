## ODEs -----
deriv(Y[, , ]) <- d_pop[i, j, k]+ d_tb[i, j, k] + d_intv[i, j, k]
dim(Y) <- c(n_tbp, n_drp, n_agp)


## Initial values -----
initial(Y[, , ]) <- Y0[i, j, k]

Y0[, , , ] <- user()
dim(Y0) <- c(n_tbp, n_drp, n_agp)


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


# TB status
I_U <- 1
I_Asym <- 2
I_Sym <- 3 # Symptomatic without care-seeking
I_ExCS <- 4 # Revisit care-seeking
I_ReCS <- 5 # Returned care-seeker after LTFU from tx
I_TxPub <- 6
I_TxPriOnPub <- 7
I_TxPriOnPri <- 8
I_TxsPub <- 9
I_TxsPriOnPub <- 10

I_FLat <- 11
I_SLat <- 12
I_RHighPub <- 13
I_RStPub <- 14
I_RHighPri <- 15
I_RStPri <- 16

## lengths -----
n_agp <- user()
n_tbp <- 16
n_drp <- 3
n_tt_demo <- length(tt_demo)
n_ttp <- length(ttp)

## birth -----

r_birth[] <- user()
dim(r_birth) <- n_tt_demo
br_t <- interpolate(tt1, r_birth, "linear")


## death + migration -----
omega[, , ] <- user()
dim(omega) <- c(n_tt1, n_agp, 2)

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


## TB -----

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