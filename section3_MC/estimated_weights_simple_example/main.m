%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example file for composite likelihood with two NK models, each with observables(x,pi,r), where
% pi is CPI-based inflation, r is the federal funds rate and x is either real per capita GDP growth or real per capita GDP detrended with an HP filter.
%The model is now estimated with multiple block MH
%and actual data from third quarter of 1954 to the end of 2015, except for the model
%with core CPI, where data starts in 1957 (2nd quarter). All parameters are common for
%all models except for the mean gdp growth rate and the parameter phi.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;





%now set up the estimation

setup_estimated_weights;

%and run the estimation

[ draws, acc_rate, log_posteriors, statedraws, individual_post_kernels] = sampling_MH( setup );



save results_equal_prior