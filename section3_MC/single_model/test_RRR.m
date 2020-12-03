clear all
clc;


sigma=1/8;
gamma=1;
rho_r=0.5;
gamma_pi=1.5;
gamma_y=0.2;
rho_a=0.8;
sigma_a=2/100;
rho_g=0.8;
sigma_g=12/100;
sigma_z=0.5/100;
sigma_lambda=4/10;
gamma_b=.2;
kappa_w=.1;
kappa_p=.1;

params=[sigma gamma rho_r gamma_pi gamma_y rho_a sigma_a rho_g sigma_g sigma_z sigma_lambda gamma_b kappa_w kappa_p]';


[A B C_final D_final add_matrices] = RRR( params,[],[] );