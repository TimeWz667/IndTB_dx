data {
  // building up PPV and Pr(Dx)
  int<lower=0> Pop[2];
  
  real tp_bac[4];
  real fp_bac[4];
  real test_tb_ssm[4];
  real test_tb_xpert[4];
  real test_nontb_ssm[4];
  real test_nontb_xpert[4];
  
  
  // Tx
  int<lower=0> N_Test_SSM_Pub[2];
  int<lower=0> N_Test_Xpert_Pub[2];
  int<lower=0> N_Test_Xpert_Eng[2];
  
  int<lower=0> N_Det_Bac[2];
  int<lower=0> N_Det_Xpert[2];
  int<lower=0> N_Det_CDx[2];
  int<lower=0> N_Det_Pub[2];
  int<lower=0> N_Det_Eng[2];
  
  int<lower=0> N_Txi_Pub[2];
  int<lower=0> N_Txi_Eng[2];
  
  real<lower=0> Drug;
  real<lower=0> Drug_Std;
  
  int<lower=0> Tx;
  int<lower=0> Tx_Pub;
  
  real<lower=0, upper=1> p_csi_pub;
  
  real<lower=0.5> dur_upper;

}
parameters {
  real<lower=0, upper=0.1> r_test[2];
  real<lower=0, upper=1> p_tb_pub;
  real<lower=0, upper=1> p_tb_pri;
  
  real<lower=0, upper=1> rp_ent;
  
  real<lower=0, upper=1> sens_cdx;
  real<lower=0, upper=1> spec_cdx;
  
  real<lower=0, upper=1> p_ava_ssm_pub[2];
  real<lower=0, upper=1> p_ava_xpert_pub[2];
  real<lower=0, upper=1> p_ava_xpert_eng[2];
  
  // real<lower=0, upper=1> p_csi_pub;
  real<lower=0.05, upper=1> p_csi_ppm;
  
  real<lower=0.5, upper=1> p_txi_pub;
  real<lower=0.5, upper=1> p_txi_eng;
  real<lower=0.5, upper=0.8> p_txi_pri;
  
  real<lower=0.04166667, upper=dur_upper> dur_pri;
  real<lower=0, upper=1> p_pri_on_pub;
  
}
transformed parameters {
  real<lower=0> r_tb_pub[2];
  real<lower=0> r_tb_eng[2];
  real<lower=0> r_tb_pri[2];
  
  real<lower=0> r_nontb_pub[2];
  real<lower=0> r_nontb_eng[2];
  real<lower=0> r_nontb_pri[2];
  
  
  real<lower=0, upper=1> alg_pub[4, 2];
  real<lower=0, upper=1> alg_eng[2, 2];
  
  real<lower=0> p_test_tb_ssm_pub[2];
  real<lower=0> p_test_tb_xpert_pub[2];
  real<lower=0> p_test_nontb_ssm_pub[2];
  real<lower=0> p_test_nontb_xpert_pub[2];
  real<lower=0, upper=1> p_tp_bac_pub[2];
  real<lower=0, upper=1> p_fp_bac_pub[2];
  real<lower=0, upper=1> p_tp_cdx_pub[2];
  real<lower=0, upper=1> p_fp_cdx_pub[2];
  
  real<lower=0> p_test_tb_xpert_eng[2];
  real<lower=0> p_test_nontb_xpert_eng[2];
  real<lower=0, upper=1> p_tp_bac_eng[2];
  real<lower=0, upper=1> p_fp_bac_eng[2];
  real<lower=0, upper=1> p_tp_cdx_eng[2];
  real<lower=0, upper=1> p_fp_cdx_eng[2];
  
  real<lower=0, upper=1> p_tp_cdx_pri[2];
  real<lower=0, upper=1> p_fp_cdx_pri[2];
  
  real<lower=0> r_test_ssm_pub[2];
  real<lower=0> r_test_xpert_pub[2];
  real<lower=0> r_test_xpert_eng[2];
  real<lower=0> r_det_bac_pub[2];
  real<lower=0> r_det_cdx_pub[2];
  
  real<lower=0> r_det_bac_eng[2];
  real<lower=0> r_det_cdx_eng[2];
  real<lower=0> r_det_cdx_pri[2];
  
  real<lower=0> r_det_pub[2];
  real<lower=0> r_det_eng[2];
  real<lower=0> r_det_bac[2];
  real<lower=0> r_det_cdx[2];

  real<lower=0> r_txi_pub[2];
  real<lower=0> r_txi_eng[2];
  real<lower=0> r_txi_pri[2];
  
  real p_tb_eng = p_tb_pri;
  
  
  real drug_time;
  real<lower=0, upper=1> p_pub;

  // Entry points
  
  for (t in 1:2) {
    r_tb_pub[t] = r_test[t] * p_tb_pub * p_csi_pub;
    r_tb_eng[t] = r_test[t] * p_tb_eng * (1 - p_csi_pub) * p_csi_ppm * rp_ent;
    r_tb_pri[t] = r_test[t] * p_tb_pri * (1 - p_csi_pub) * (1 - p_csi_ppm) * rp_ent;
    
    r_nontb_pub[t] = r_test[t] * (1 - p_tb_pub) * p_csi_pub;
    r_nontb_eng[t] = r_test[t] * (1 - p_tb_eng) * (1 - p_csi_pub) * p_csi_ppm * rp_ent;
    r_nontb_pri[t] = r_test[t] * (1 - p_tb_pri) * (1 - p_csi_pub) * (1 - p_csi_ppm) * rp_ent;
    
    alg_pub[1, t] =      p_ava_ssm_pub[t]  *      p_ava_xpert_pub[t];
    alg_pub[2, t] =      p_ava_ssm_pub[t]  * (1 - p_ava_xpert_pub[t]);
    alg_pub[3, t] = (1 - p_ava_ssm_pub[t]) *      p_ava_xpert_pub[t];
    alg_pub[4, t] = (1 - p_ava_ssm_pub[t]) * (1 - p_ava_xpert_pub[t]);
    
    alg_eng[1, t] = p_ava_xpert_eng[t];
    alg_eng[2, t] = (1 - p_ava_xpert_eng[t]);
    
    // Public
    p_test_tb_ssm_pub[t] = 0;
    p_test_tb_xpert_pub[t] = 0;
    p_test_nontb_ssm_pub[t] = 0;
    p_test_nontb_xpert_pub[t] = 0;
    p_tp_bac_pub[t] = 0;
    p_fp_bac_pub[t] = 0;
    p_tp_cdx_pub[t] = 0;
    p_fp_cdx_pub[t] = 0;
      
    // Private, engaged
    p_test_tb_xpert_eng[t] = 0;
    p_test_nontb_xpert_eng[t] = 0;
    p_tp_bac_eng[t] = 0;
    p_fp_bac_eng[t] = 0;
    p_tp_cdx_eng[t] = 0;
    p_fp_cdx_eng[t] = 0;
    
    for(i in 1:4) {
      p_test_tb_ssm_pub[t] +=      alg_pub[i, t] * test_tb_ssm[i];
      p_test_tb_xpert_pub[t] +=    alg_pub[i, t] * test_tb_xpert[i];
      p_test_nontb_ssm_pub[t] +=   alg_pub[i, t] * test_nontb_ssm[i];
      p_test_nontb_xpert_pub[t] += alg_pub[i, t] * test_nontb_xpert[i];
  
      p_tp_bac_pub[t] += alg_pub[i, t] * tp_bac[i];
      p_fp_bac_pub[t] += alg_pub[i, t] * fp_bac[i];
      
      p_tp_cdx_pub[t] += alg_pub[i, t] * (1 - tp_bac[i]) * sens_cdx;
      p_fp_cdx_pub[t] += alg_pub[i, t] * (1 - fp_bac[i]) * (1 - spec_cdx);
      
      if (i >= 3) {
        p_test_tb_xpert_eng[t] +=    alg_eng[i-2, t] * test_tb_xpert[i];
        p_test_nontb_xpert_eng[t] += alg_eng[i-2, t] * test_nontb_xpert[i];
        
        p_tp_bac_eng[t] += alg_eng[i-2, t] * tp_bac[i];
        p_fp_bac_eng[t] += alg_eng[i-2, t] * fp_bac[i];
        
        p_tp_cdx_eng[t] += alg_eng[i-2, t] * (1 - tp_bac[i]) * sens_cdx;
        p_fp_cdx_eng[t] += alg_eng[i-2, t] * (1 - fp_bac[i]) * (1 - spec_cdx);
      }
    }
  
    p_tp_cdx_pri[t] = sens_cdx;
    p_fp_cdx_pri[t] = (1 - spec_cdx);
    
  
    // Map to notification data
    r_test_ssm_pub[t] =   r_tb_pub[t] * p_test_tb_ssm_pub[t]   + r_nontb_pub[t] * p_test_nontb_ssm_pub[t];
    r_test_xpert_pub[t] = r_tb_pub[t] * p_test_tb_xpert_pub[t] + r_nontb_pub[t] * p_test_nontb_xpert_pub[t];
    r_det_bac_pub[t] =    r_tb_pub[t] * p_tp_bac_pub[t]        + r_nontb_pub[t] * p_fp_bac_pub[t];
    r_det_cdx_pub[t] =    r_tb_pub[t] * p_tp_cdx_pub[t]        + r_nontb_pub[t] * p_fp_cdx_pub[t];
    
    r_test_xpert_eng[t] = r_tb_eng[t] * p_test_tb_xpert_eng[t] + r_nontb_eng[t] * p_test_nontb_xpert_eng[t];
    r_det_bac_eng[t] =    r_tb_eng[t] * p_tp_bac_eng[t]        + r_nontb_eng[t] * p_fp_bac_eng[t];
    r_det_cdx_eng[t] =    r_tb_eng[t] * p_tp_cdx_eng[t]        + r_nontb_eng[t] * p_fp_cdx_eng[t];
    
    r_det_pub[t] = r_det_bac_pub[t] + r_det_cdx_pub[t];
    r_det_eng[t] = r_det_bac_eng[t] + r_det_cdx_eng[t];
    r_det_bac[t] = r_det_bac_pub[t] + r_det_bac_eng[t];
    r_det_cdx[t] = r_det_cdx_pub[t] + r_det_cdx_eng[t];
    
    r_det_cdx_pri[t] = r_tb_pri[t] * p_tp_cdx_pri[t] + r_nontb_pri[t] * p_fp_cdx_pri[t];
    
    r_txi_pub[t] = r_det_pub[t] * p_txi_pub;
    r_txi_eng[t] = r_det_eng[t] * p_txi_eng;
    r_txi_pri[t] = r_det_cdx_pri[t] * p_txi_pri;

    
  }

  drug_time = (r_txi_pri[1] + r_txi_eng[1] * (1 - p_pri_on_pub)) * dur_pri;
  
  p_pub = sum(r_txi_pub) / (sum(r_txi_pub) + sum(r_txi_eng) + sum(r_txi_pri));
  
}
model {
  p_tb_pub ~ uniform(0, 1);
  p_tb_pri ~ uniform(0, 1);
  p_csi_ppm ~ uniform(0, 1);
  p_pri_on_pub ~ beta(1.5, 3.5);
  p_txi_pri ~ uniform(0.5, 0.8);
  
  for (i in 1:2) {
    target += binomial_lpmf(N_Test_SSM_Pub[i] | Pop, r_test_ssm_pub[i]);
    target += binomial_lpmf(N_Test_Xpert_Pub[i] | Pop, r_test_xpert_pub[i]);
    target += binomial_lpmf(N_Test_Xpert_Eng[i] | Pop, r_test_xpert_eng[i]);
    target += binomial_lpmf(N_Txi_Pub[i] | N_Det_Pub[i], p_txi_pub);
    target += binomial_lpmf(N_Txi_Eng[i] | N_Det_Eng[i], p_txi_eng);
    
    target += binomial_lpmf(N_Det_Pub[i] | Pop[i], r_det_pub[i]);
    target += binomial_lpmf(N_Det_Eng[i] | Pop[i], r_det_eng[i]);
    
    target += binomial_lpmf(N_Det_Bac[i] | Pop[i], r_det_bac[i]);
    target += binomial_lpmf(N_Det_CDx[i] | Pop[i], r_det_cdx[i]);

  }

  target += binomial_lpmf(Tx_Pub | Tx, p_pub);
  target += normal_lpdf(Drug | drug_time, Drug_Std);

}
generated quantities {
  real<lower=0, upper=1> ppv_pub[2];
  real<lower=0, upper=1> ppv_eng[2];
  real<lower=0, upper=1> ppv_pri[2];
  
  real<lower=0, upper=1> p_under[2];
  real<lower=0> tp_pri_drug[2];
  real<lower=0> tp_pri_drug_time[2];
  real<lower=0> tp_pri_txi[2];
    
  
  for (i in 1:2) {
    ppv_pub[i] = r_tb_pub[i] * (p_tp_bac_pub[i] + p_tp_cdx_pub[i]) / r_det_pub[i];
    ppv_eng[i] = r_tb_eng[i] * (p_tp_bac_eng[i] + p_tp_cdx_eng[i]) / r_det_eng[i];
    ppv_pri[i] = r_tb_pri[i] * p_tp_cdx_pri[i] / r_det_cdx_pri[i];
    
    
    tp_pri_drug[i]  = Pop[i] * r_tb_eng[i] * (p_tp_bac_eng[i] + p_tp_cdx_eng[i]) * (1 - p_pri_on_pub) * p_txi_eng;
    tp_pri_drug[i] += Pop[i] * r_tb_pri[i] * p_tp_cdx_pri[i] * p_txi_pri;
    
    tp_pri_txi[i]  = Pop[i] * r_tb_eng[i] * (p_tp_bac_eng[i] + p_tp_cdx_eng[i]) * p_txi_eng;
    tp_pri_txi[i] += Pop[i] * r_tb_pri[i] * p_tp_cdx_pri[i] * p_txi_pri;
    
    tp_pri_drug_time[i] = tp_pri_drug[i] * dur_pri;
    p_under[i] = r_det_cdx_pri[i] / (r_det_pub[i] + r_det_eng[i] + r_det_cdx_pri[i]);
    
  }

}
