%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example file for composite likelihood with four NK models, each with observables(delta log x,pi,r), where
% pi is either CPI-based inflation, core CPI based inflation, PCE-based inflation or 
% GDP-deflator-based inflation. The model is now estimated with multiple block MH
%and actual data from third quarter of 1954 to the end of 2015, except for the model
%with core CPI, where data starts in 1957 (2nd quarter). All parameters are common for
%all models except for mean inflation rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;





%now set up the estimation

setup_DSGE_data;

%and run the estimation

[ draws, acc_rate, log_posteriors, statedraws, individual_post_kernels] = sampling_MH( setup );


% figure;
% subplot(4,5,1)
% hist(draws(1,:),50)
% title('\sigma')
% subplot(4,5,2)
% hist(draws(2,:),50)
% title('\gamma')
% subplot(4,5,3)
% hist(draws(3,:),50)
% title('\rho_r')
% subplot(4,5,4)
% hist(draws(4,:),50)
% title('\gamma_{\pi}')
% subplot(4,5,5)
% hist(draws(5,:),50)
% title('\gamma_y')
% subplot(4,5,6)
% hist(draws(6,:),50)
% title('\rho_a')
% subplot(4,5,7)
% hist(draws(7,:),50)
% title('\sigma_a')
% subplot(4,5,8)
% hist(draws(8,:),50)
% title('\rho_g')
% subplot(4,5,9)
% hist(draws(9,:),50)
% title('\sigma_g')
% subplot(4,5,10)
% hist(draws(10,:),50)
% title('\sigma_z')
% subplot(4,5,11)
% hist(draws(11,:),50)
% title('\sigma_{\lambda}')
% subplot(4,5,12)
% hist(draws(11,:),50)
% title('\gamma_b')
% subplot(4,5,13)
% hist(draws(11,:),50)
% title('\kappa_w')
% subplot(4,5,13)
% hist(draws(11,:),50)
% title('\kappa_p')

save results2
