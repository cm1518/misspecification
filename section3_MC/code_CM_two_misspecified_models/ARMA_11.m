function [A B C D add_matrices] = ARMA_11( params,setup,data )
%state space matrices for ARMA_11 model
sigma=params(1);
theta=params(2);
rho=params(3);


add_matrices=[];

A=[1 0 0];
B=0;
C=[rho theta 0;0 0 0;0 1 0];
D=[sigma^2 0 0;sigma^2 sigma^2 0;0 0 0];