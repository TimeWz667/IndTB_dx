data {
  // building up PPV and Pr(Dx)
  int<lower=0> Pop;
  
  // ITR
  int<lower=0> N_Pre_SSM;
  int<lower=0> N_Pre_Xpert;
  
  int<lower=0> N_Det_SSM_Pub;
  int<lower=0> N_Det_Xpert_Pub;
  int<lower=0> N_Det_Xpert_Eng;
  
  int<lower=0> N_Det_CDX_Pub;
  int<lower=0> N_Det_CDX_Eng;
}
parameters {
  real<lower=0, upper=0.1> r_test;
  real<lower=0, upper=1> p_tb;
  
  real<lower=0, upper=1> p_pub;
  real<lower=0, upper=1> p_path_ssm;
  real<lower=0, upper=1> p_xpert_sn;
  real<lower=0.05, upper=0.5> p_xpert_eng;
  real<lower=0.15, upper=0.5> p_scanty;

  real<lower=0.63, upper=0.65> sens_s;
  real<lower=0.82, upper=0.88> sens_x;
  real<lower=0.63, upper=0.65> sens_x_sn;
  real<lower=0.97, upper=0.99> spec_s;
  real<lower=0.97, upper=0.99> spec_x;
  real<lower=0.97, upper=0.99> spec_x_sn;
  
  real<lower=0.95, upper=1> sens_cdx;
  real<lower=0.72, upper=0.79> spec_cdx;
  
  real<lower=0, upper=1> sens_cdx_bn;
  real<lower=0.5, upper=1> spec_cdx_bn;
}
transformed parameters {
  real r_access_s_pub = r_test * p_pub * p_path_ssm * (1 - p_xpert_sn);
  real r_access_x_pub = r_test * p_pub * (1 - p_path_ssm);
  real r_access_sx_pub = r_test * p_pub * p_path_ssm * p_xpert_sn;
  
  real r_access_x_eng = r_test * (1 - p_pub) * p_xpert_eng;
  real r_access_cdx_eng = r_test * (1 - p_pub) * (1 - p_xpert_eng);
  
  // public system
  real r_tpb_s_pub = r_access_s_pub *      p_tb  *  (1 - p_scanty) *      sens_s;
  real r_fpb_s_pub = r_access_s_pub * (1 - p_tb) *  (1 - p_scanty) * (1 - spec_s);
  
  real r_tpb_x_pub = r_access_x_pub *      p_tb  *  (1 - p_scanty) *      sens_x;
  real r_fpb_x_pub = r_access_x_pub * (1 - p_tb) *  (1 - p_scanty) * (1 - spec_x);
  
  real r_tps_sx_pub = r_access_sx_pub *      p_tb  *  (1 - p_scanty) *      sens_s;
  real r_fps_sx_pub = r_access_sx_pub * (1 - p_tb) *  (1 - p_scanty) * (1 - spec_s);
  real r_tpx_sx_pub = r_access_sx_pub *      p_tb  *  (1 - p_scanty) * (1 - sens_s) *      sens_x_sn;
  real r_fpx_sx_pub = r_access_sx_pub * (1 - p_tb) *  (1 - p_scanty) *      spec_s  * (1 - spec_x_sn);
  real r_tpb_sx_pub = r_tps_sx_pub + r_tpx_sx_pub;
  real r_fpb_sx_pub = r_fps_sx_pub + r_fpx_sx_pub;
  
  real r_tpb_pub = r_tpb_s_pub + r_tpb_x_pub + r_tpb_sx_pub;
  real r_fpb_pub = r_fpb_s_pub + r_fpb_x_pub + r_fpb_sx_pub;
  real r_tpc_pub = ((r_access_s_pub + r_access_x_pub + r_access_sx_pub) *      p_tb  - r_tpb_pub) *      sens_cdx_bn;
  real r_fpc_pub = ((r_access_s_pub + r_access_x_pub + r_access_sx_pub) * (1 - p_tb) - r_fpb_pub) * (1 - spec_cdx_bn);

  // private system
  real r_tpb_x_eng = r_access_x_eng *      p_tb  *  (1 - p_scanty) *      sens_x;
  real r_fpb_x_eng = r_access_x_eng * (1 - p_tb) *  (1 - p_scanty) * (1 - spec_x);
  
  real r_tpc_eng = ((r_access_x_eng) *      p_tb  - r_tpb_x_eng) *      sens_cdx_bn;
  real r_fpc_eng = ((r_access_x_eng) * (1 - p_tb) - r_fpb_x_eng) * (1 - spec_cdx_bn);
  
  r_tpc_eng += r_access_cdx_eng *      p_tb  *      sens_cdx;
  r_fpc_eng += r_access_cdx_eng * (1 - p_tb) * (1 - spec_cdx);
  
}
model {
  p_tb ~ uniform(0, 1);
  p_xpert_eng ~ uniform(0.05, 0.5);
  
  target += binomial_lpmf(N_Pre_SSM | Pop, r_access_s_pub + r_access_sx_pub);
  target += binomial_lpmf(N_Pre_Xpert | Pop, r_access_x_pub + r_access_x_eng);
  
  target += binomial_lpmf(N_Det_SSM_Pub | Pop, r_tpb_s_pub + r_fpb_s_pub + r_tps_sx_pub + r_fps_sx_pub);
  target += binomial_lpmf(N_Det_Xpert_Pub | Pop, r_tpb_x_pub + r_fpb_x_pub + r_tpx_sx_pub + r_fpx_sx_pub);
  target += binomial_lpmf(N_Det_Xpert_Eng | Pop, r_tpb_x_eng + r_fpb_x_eng);
  
  target += binomial_lpmf(N_Det_CDX_Pub | Pop, r_tpc_pub + r_fpc_pub);
  target += binomial_lpmf(N_Det_CDX_Eng | Pop, r_tpc_eng + r_fpc_eng);

}
generated quantities {
  real ppv_pub = (r_tpb_pub + r_tpc_pub) / (r_tpb_pub + r_tpc_pub + r_fpb_pub + r_fpc_pub);
  real ppv_eng = (r_tpb_x_eng + r_tpc_eng) / (r_tpb_x_eng + r_tpc_eng + r_fpb_x_eng + r_fpc_eng);
  real ppv_pri = p_tb * sens_cdx / (p_tb * sens_cdx + (1 - p_tb) * (1 - spec_cdx));
  
  real pdx_pub = (r_tpb_pub + r_tpc_pub) / (r_access_s_pub + r_access_x_pub + r_access_sx_pub) / p_tb;
  real pdx_eng = (r_tpb_x_eng + r_tpc_eng) / (r_access_x_eng + r_access_cdx_eng) / p_tb;
  real pdx_pri = sens_cdx;
}
