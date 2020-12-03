function [A B C D add_matrices] = stocks( params,setup,data )
%state space matrices for MA_2 model
sigma=params(1);
theta=params(2);
beta=params(3);

add_matrices=[];

A=[1-beta*theta -theta];
B=0;
C=[0 0;1 0];
D=[sigma^2 0;0 0];