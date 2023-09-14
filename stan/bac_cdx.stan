data {
  real tp_bac[4];
  real fp_bac[4];
  real test_tb_ssm[4];
  real test_tb_xpert[4];
  real test_nontb_ssm[4];
  real test_nontb_xpert[4];
  

  int<lower=0> Pop;
  int<lower=0> N_Test_SSM;
  int<lower=0> N_Test_Xpert;
  
  int<lower=0> N_Det_Bac;
  int<lower=0> N_Det_Xpert;
  int<lower=0> N_Det_CDx;

}
parameters {
  real<lower=0, upper=0.1> r_test;
  real<lower=0, upper=1> p_tb;
  
  real<lower=0, upper=1> sens_cdx;
  real<lower=0, upper=1> spec_cdx;
  
  real<lower=0, upper=1> p_ava_ssm;
  real<lower=0, upper=1> p_ava_xpert;
  real<lower=0, upper=1> p_ava_xpert_ssm;
}
transformed parameters {
  real<lower=0> r_tb;
  real<lower=0> r_nontb;

  real<lower=0, upper=1> alg[4];

  
  real<lower=0> p_test_tb_ssm = 0;
  real<lower=0> p_test_tb_xpert = 0;
  real<lower=0> p_test_nontb_ssm = 0;
  real<lower=0> p_test_nontb_xpert = 0;

  real<lower=0, upper=1> p_tp_bac = 0;
  real<lower=0, upper=1> p_tp_cdx = 0;
  real<lower=0, upper=1> p_fp_bac = 0;
  real<lower=0, upper=1> p_fp_cdx = 0;
  
  real<lower=0> r_test_ssm = 0;
  real<lower=0> r_test_xpert = 0;

  real<lower=0> r_det_bac = 0;
  real<lower=0> r_det_cdx = 0;

  // Entry points
  
  r_tb = r_test * p_tb;
  r_nontb = r_test * (1 - p_tb);
  
  alg[1] =      p_ava_ssm  *      p_ava_xpert_ssm;
  alg[2] =      p_ava_ssm  * (1 - p_ava_xpert_ssm);
  alg[3] = (1 - p_ava_ssm) *      p_ava_xpert;
  alg[4] = (1 - p_ava_ssm) * (1 - p_ava_xpert);
  
  for(i in 1:4) {
    p_test_tb_ssm +=      alg[i] * test_tb_ssm[i];
    p_test_tb_xpert +=    alg[i] * test_tb_xpert[i];
    p_test_nontb_ssm +=   alg[i] * test_nontb_ssm[i];
    p_test_nontb_xpert += alg[i] * test_nontb_xpert[i];

    p_tp_bac += alg[i] * tp_bac[i];
    p_fp_bac += alg[i] * fp_bac[i];
    
    p_tp_cdx += alg[i] * (1 - tp_bac[i]) * sens_cdx;
    p_fp_cdx += alg[i] * (1 - fp_bac[i]) * (1 - spec_cdx);
  }
  

  r_test_ssm =   r_tb * p_test_tb_ssm   + r_nontb * p_test_nontb_ssm;
  r_test_xpert = r_tb * p_test_tb_xpert + r_nontb * p_test_nontb_xpert;
  r_det_bac =    r_tb * p_tp_bac        + r_nontb * p_fp_bac;
  r_det_cdx =    r_tb * p_tp_cdx        + r_nontb * p_fp_cdx;
}
model {
  p_tb ~ uniform(0, 1);
  
  target += binomial_lpmf(N_Test_SSM | Pop, r_test_ssm);
  target += binomial_lpmf(N_Test_Xpert | Pop, r_test_xpert);
  
  target += binomial_lpmf(N_Det_Bac | Pop, r_det_bac);
  target += binomial_lpmf(N_Det_CDx | Pop, r_det_cdx);

}
