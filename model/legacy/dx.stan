data {
  // building up PPV and Pr(Dx)
  int<lower=0> Pop;
  
  real<lower=0, upper=1> Sens_SSM;
  real<lower=0, upper=1> Sens_Xpert;
  real<lower=0, upper=1> Sens_XpertSN;
  
  real<lower=0, upper=1> Spec_SSM;
  real<lower=0, upper=1> Spec_Xpert;
  real<lower=0, upper=1> Spec_XpertSN;
  
  
  // ITR
  int<lower=0> N_Pre_SSM;
  int<lower=0> N_Pre_Xpert;
  
  int<lower=0> N_Det_SSM;
  int<lower=0> N_Det_Xpert;
  
  int<lower=0> N_Det_CDX;
  
  int<lower=0> N_Det_Pub;
  int<lower=0> N_Det_Eng;
  
  real<lower=0, upper=1> PrEng;
}
parameters {
  real<lower=0, upper=0.1> r_test;
  real<lower=0, upper=1> p_tb;
  real<lower=0, upper=1> sens_cdx;
  real<lower=0, upper=1> spec_cdx;
  
  real<lower=0, upper=1> p_test_bac;
  real<lower=0, upper=1> p_test_ssm_bac;
  real<lower=0, upper=1> p_xpert_sn;
  real<lower=0.15, upper=0.5> p_scanty;

}
transformed parameters {
  real r_ent_ssm = r_test * p_test_ssm_bac;
  real r_ent_xpert = r_test * (1- p_test_ssm_bac);
  real r_ent_cdx = r_test * (1 - p_test_bac);
  
  real r_tp_ssm = r_ent_ssm *      p_tb  *  (1 - p_scanty) *      Sens_SSM;
  real r_fp_ssm = r_ent_ssm * (1 - p_tb) *  (1 - p_scanty) * (1 - Spec_SSM);
  
  real r_tp_xpert_sn = r_ent_ssm *      p_tb  *  (1 - p_scanty) * (1 - Sens_SSM) * p_xpert_sn *      Sens_XpertSN;
  real r_fp_xpert_sn = r_ent_ssm * (1 - p_tb) *  (1 - p_scanty) *      Spec_SSM  * p_xpert_sn * (1 - Spec_XpertSN);
  
  real r_tp_xpert = r_ent_xpert *      p_tb  *  (1 - p_scanty) *      Sens_Xpert;
  real r_fp_xpert = r_ent_xpert * (1 - p_tb) *  (1 - p_scanty) * (1 - Spec_Xpert);

  real r_tp_cdx = ((r_ent_ssm + r_ent_xpert) *      p_tb  - r_tp_ssm - r_tp_xpert_sn - r_tp_xpert) *      sens_cdx;
  real r_fp_cdx = ((r_ent_ssm + r_ent_xpert) * (1 - p_tb) - r_fp_ssm - r_fp_xpert_sn - r_fp_xpert) * (1 - spec_cdx);
  
  real r_det_ssm = r_tp_ssm + r_fp_ssm;
  real r_det_xpert = r_tp_xpert + r_fp_xpert + r_tp_xpert_sn + r_fp_xpert_sn;
  real r_det_cdx = r_tp_cdx + r_fp_cdx;
}
model {
  p_tb ~ uniform(0, 1);
  
  target += binomial_lpmf(N_Pre_SSM | Pop, r_ent_ssm);
  target += binomial_lpmf(N_Pre_Xpert | Pop, r_ent_xpert);
  
  target += binomial_lpmf(N_Det_SSM | Pop, r_det_ssm);
  target += binomial_lpmf(N_Det_Xpert | Pop, r_det_xpert);
  target += binomial_lpmf(N_Det_CDX | Pop, r_det_cdx);

}
generated quantities {
  real<lower=0, upper=1> ppv_ssm = r_tp_ssm / r_det_ssm;
  real<lower=0, upper=1> ppv_xpert = (r_tp_xpert + r_tp_xpert_sn) / r_det_xpert;
  real<lower=0, upper=1> ppv_cdx = r_tp_cdx / r_det_cdx;
  

  
}
