%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for technology shock application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;





%set up the estimation

setup_technology_QE;

%and run the estimation

[ draws, acc_rate, log_posteriors, statedraws, individual_post_kernels] = sampling_MH( setup );



save results_technology_QE

var_decomp_plot
