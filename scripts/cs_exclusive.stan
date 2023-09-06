data {
  // building up PPV and Pr(Dx)
  int<lower=0> Pop;
  
  real tp_bac[4];
  real fp_bac[4];
  real test_tb_ssm[4];
  real test_tb_naat[4];
  real test_nontb_ssm[4];
  real test_nontb_naat[4];
  
  
  // Tx
  int<lower=0> N_Test_SSM_Pub;
  int<lower=0> N_Test_NAAT_Pub;
  int<lower=0> N_Test_NAAT_Eng;
  
  int<lower=0> N_DetBac;
  int<lower=0> N_DetCDx;
  int<lower=0> N_Det_Pub;
  int<lower=0> N_Det_Eng;
  // int<lower=0> N_DetBac_Pub;
  // int<lower=0> N_DetCDx_Pub;
  // int<lower=0> N_DetBac_Eng;
  // int<lower=0> N_DetCDx_Eng;
  
  int<lower=0> N_Txi_Pub;
  int<lower=0> N_Txi_Eng;
  
  real<lower=0> Drug;
  real<lower=0> Drug_Std;
  
  int<lower=0> Tx;
  int<lower=0> Tx_Pub;
  
  real<lower=0, upper=1> p_csi_pub;
  
  real<lower=0.5> dur_upper;

}
parameters {
  real<lower=0, upper=0.1> r_test;
  real<lower=0, upper=1> p_tb;
  real<lower=0, upper=1> sens_cdx;
  real<lower=0, upper=1> spec_cdx;
  
  real<lower=0, upper=1> p_ava_ssm_pub;
  real<lower=0, upper=1> p_ava_naat_pub;
  real<lower=0, upper=1> p_ava_naat_eng;
  
  // real<lower=0, upper=1> p_csi_pub;
  real<lower=0.05, upper=1> p_csi_ppm;
  
  real<lower=0.5, upper=1> p_txi_pub;
  real<lower=0.5, upper=1> p_txi_eng;
  real<lower=0.5, upper=0.8> p_txi_pri;
  
  real<lower=0.04166667, upper=dur_upper> dur_pri;
  real<lower=0, upper=1> p_pri_on_pub;
  
}
transformed parameters {
  real<lower=0, upper=1> alg_pub[4];
  real<lower=0, upper=1> alg_eng[2];
  
  real<lower=0> p_test_tb_ssm_pub;
  real<lower=0> p_test_tb_naat_pub;
  real<lower=0> p_test_nontb_ssm_pub;
  real<lower=0> p_test_nontb_naat_pub;
  real<lower=0, upper=1> p_tp_bac_pub;
  real<lower=0, upper=1> p_fp_bac_pub;
  real<lower=0, upper=1> p_tp_cdx_pub;
  real<lower=0, upper=1> p_fp_cdx_pub;
  
  real<lower=0> p_test_tb_naat_eng;
  real<lower=0> p_test_nontb_naat_eng;
  real<lower=0, upper=1> p_tp_bac_eng;
  real<lower=0, upper=1> p_fp_bac_eng;
  real<lower=0, upper=1> p_tp_cdx_eng;
  real<lower=0, upper=1> p_fp_cdx_eng;
  
  real<lower=0, upper=1> p_tp_cdx_pri;
  real<lower=0, upper=1> p_fp_cdx_pri;
  
  real<lower=0> r_test_ssm_pub;
  real<lower=0> r_test_naat_pub;
  real<lower=0> r_test_naat_eng;
  real<lower=0> r_det_bac_pub;
  real<lower=0> r_det_cdx_pub;
  real<lower=0> r_det_bac_eng;
  real<lower=0> r_det_cdx_eng;
  real<lower=0> r_det_cdx_pri;

  real<lower=0> r_txi_pub;
  real<lower=0> r_txi_eng;
  real<lower=0> r_txi_pri;
  
  
  real drug_time;
  real<lower=0, upper=1> p_pub;

  alg_pub[1] = 0;
  alg_pub[2] = p_ava_ssm_pub;
  alg_pub[3] = (1 - p_ava_ssm_pub) * p_ava_naat_pub;
  alg_pub[4] = (1 - p_ava_ssm_pub) * (1 - p_ava_naat_pub);
  
  alg_eng[1] = p_ava_naat_eng;
  alg_eng[2] = (1 - p_ava_naat_eng);


  // Public
  p_test_tb_ssm_pub = 0;
  p_test_tb_naat_pub = 0;
  p_test_nontb_ssm_pub = 0;
  p_test_nontb_naat_pub = 0;
  p_tp_bac_pub = 0;
  p_fp_bac_pub = 0;
  p_tp_cdx_pub = 0;
  p_fp_cdx_pub = 0;
  
  for(i in 1:4) {
    p_test_tb_ssm_pub +=     p_tb       * alg_pub[i] * test_tb_ssm[i];
    p_test_tb_naat_pub +=    p_tb       * alg_pub[i] * test_tb_naat[i];
    p_test_nontb_ssm_pub +=  (1 - p_tb) * alg_pub[i] * test_nontb_ssm[i];
    p_test_nontb_naat_pub += (1 - p_tb) * alg_pub[i] * test_nontb_naat[i];

    p_tp_bac_pub +=      p_tb  * alg_pub[i] * tp_bac[i];
    p_fp_bac_pub += (1 - p_tb) * alg_pub[i] * fp_bac[i];
    
    p_tp_cdx_pub +=      p_tb  * alg_pub[i] * (1 - tp_bac[i]) * sens_cdx;
    p_fp_cdx_pub += (1 - p_tb) * alg_pub[i] * (1 - fp_bac[i]) * (1 - spec_cdx);
  }
  
  
  
  // Private, engaged
  p_test_tb_naat_eng = 0;
  p_test_nontb_naat_eng = 0;
  p_tp_bac_eng = 0;
  p_fp_bac_eng = 0;
  p_tp_cdx_eng = 0;
  p_fp_cdx_eng = 0;
  
  for(i in 3:4) {
    p_test_tb_naat_eng +=    p_tb       * alg_eng[i-2] * test_tb_naat[i];
    p_test_nontb_naat_eng += (1 - p_tb) * alg_eng[i-2] * test_nontb_naat[i];

    p_tp_bac_eng +=      p_tb  * alg_eng[i-2] * tp_bac[i];
    p_fp_bac_eng += (1 - p_tb) * alg_eng[i-2] * fp_bac[i];
    
    p_tp_cdx_eng +=      p_tb  * alg_eng[i-2] * (1 - tp_bac[i]) * sens_cdx;
    p_fp_cdx_eng += (1 - p_tb) * alg_eng[i-2] * (1 - fp_bac[i]) * (1 - spec_cdx);
  }

  
  
  // Private, unengaged
  p_tp_cdx_pri =      p_tb * sens_cdx;
  p_fp_cdx_pri = (1 - p_tb) * (1 - spec_cdx);
  

  // Map to notification data
  r_test_ssm_pub = r_test * p_csi_pub * (p_test_tb_ssm_pub + p_test_nontb_ssm_pub);
  r_test_naat_pub = r_test * p_csi_pub * (p_test_tb_naat_pub + p_test_nontb_naat_pub);
  r_test_naat_eng = r_test * (1 - p_csi_pub) * p_csi_ppm * (p_test_tb_naat_eng + p_test_nontb_naat_eng);
  r_det_bac_pub = r_test * p_csi_pub * (p_tp_bac_pub + p_fp_bac_pub);
  r_det_cdx_pub = r_test * p_csi_pub * (p_tp_cdx_pub + p_fp_cdx_pub);
  r_det_bac_eng = r_test * (1 - p_csi_pub) * p_csi_ppm * (p_tp_bac_eng + p_fp_bac_eng);
  r_det_cdx_eng = r_test * (1 - p_csi_pub) * p_csi_ppm * (p_tp_cdx_eng + p_fp_cdx_eng);
  r_det_cdx_pri = r_test * (1 - p_csi_pub) * (1 - p_csi_ppm) * (p_tp_cdx_pri + p_fp_cdx_pri);

  r_txi_pub = (r_det_bac_pub + r_det_cdx_pub) * p_txi_pub;
  r_txi_eng = (r_det_bac_eng + r_det_cdx_eng) * p_txi_eng;
  r_txi_pri = r_det_cdx_pri * p_txi_pri;


  drug_time = (r_txi_pri + r_txi_eng * (1 - p_pri_on_pub)) * dur_pri;
  
  p_pub = r_txi_pub / (r_txi_pub + r_txi_eng + r_txi_pri);
  
}
model {
  p_tb ~ uniform(0, 1);
  
  target += binomial_lpmf(N_Test_SSM_Pub | Pop, r_test_ssm_pub);
  target += binomial_lpmf(N_Test_NAAT_Pub | Pop, r_test_naat_pub);
  target += binomial_lpmf(N_Test_NAAT_Eng | Pop, r_test_naat_eng);
  // target += binomial_lpmf(N_DetBac_Pub | Pop, r_det_bac_pub);
  // target += binomial_lpmf(N_DetCDx_Pub | Pop, r_det_cdx_pub);
  // target += binomial_lpmf(N_DetBac_Eng | Pop, r_det_bac_eng);
  // target += binomial_lpmf(N_DetCDx_Eng | Pop, r_det_cdx_eng);
  target += binomial_lpmf(N_DetBac | Pop, r_det_bac_pub + r_det_bac_eng);
  target += binomial_lpmf(N_DetCDx | Pop, r_det_cdx_pub + r_det_cdx_eng);
  target += binomial_lpmf(N_Det_Pub | Pop, r_det_bac_pub + r_det_cdx_pub);
  target += binomial_lpmf(N_Det_Eng | Pop, r_det_bac_eng + r_det_cdx_eng);
  
  p_csi_pub ~ uniform(0, 1);
  p_csi_ppm ~ uniform(0, 1);
  p_pri_on_pub ~ beta(1.5, 3.5);
  p_txi_pri ~ uniform(0.5, 0.8);
    
  target += binomial_lpmf(N_Txi_Pub | N_Det_Pub, p_txi_pub);
  target += binomial_lpmf(N_Txi_Eng | N_Det_Eng, p_txi_eng);
    
  target += binomial_lpmf(Tx_Pub | Tx, p_pub);
  
  target += normal_lpdf(Drug | drug_time, Drug_Std);

}
generated quantities {
  real<lower=0, upper=1> ppv_pub;
  real<lower=0, upper=1> ppv_eng;
  real<lower=0, upper=1> ppv_pri;
  
  real<lower=0, upper=1> p_under;
  real<lower=0> tp_pri_drug;
  real<lower=0> tp_pri_drug_time;
  real<lower=0> tp_pri_txi;
  
  ppv_pub = (p_tp_bac_pub + p_tp_cdx_pub) / (p_tp_bac_pub + p_tp_cdx_pub + p_fp_bac_pub + p_fp_cdx_pub);
  ppv_eng = (p_tp_bac_eng + p_tp_cdx_eng) / (p_tp_bac_eng + p_tp_cdx_eng + p_fp_bac_eng + p_fp_cdx_eng);
  ppv_pri = (p_tp_cdx_pri) / (p_tp_cdx_pri + p_fp_cdx_pri);
  
  tp_pri_drug = Pop * r_test * (1 - p_csi_pub) * ((1 - p_csi_ppm) * p_tp_cdx_pri * p_txi_pri + p_csi_ppm * (p_tp_bac_eng + p_tp_cdx_eng) * p_txi_eng * (1 - p_pri_on_pub));
  tp_pri_txi = Pop * r_test * (1 - p_csi_pub) * ((1 - p_csi_ppm) * p_tp_cdx_pri * p_txi_pri + p_csi_ppm * (p_tp_bac_eng + p_tp_cdx_eng) * p_txi_eng);
  tp_pri_drug_time = tp_pri_drug * dur_pri;
  
  p_under = r_det_cdx_pri / (r_det_bac_pub + r_det_cdx_pub + r_det_bac_eng + r_det_cdx_eng + r_det_cdx_pri);
  
}
