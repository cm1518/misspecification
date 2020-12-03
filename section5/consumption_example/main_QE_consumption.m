%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code estimates the consumption example with 5 models.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;





%set up the estimation

setup_consumption_QE;

%and run the estimation

[ draws, acc_rate, log_posteriors, statedraws, individual_post_kernels] = sampling_MH( setup );



save results_QE_consumption

plot_model_weights_5_models %plotting the prior and posterior model weights