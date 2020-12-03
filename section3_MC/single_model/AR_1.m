function [A B C D add_matrices] = AR_1( params,setup,data )
%state space matrices for AR(1) model
sigma=params(1);
rho=params(2);

add_matrices=[];

A=1;
B=0;
C=rho;
D=sigma^2;