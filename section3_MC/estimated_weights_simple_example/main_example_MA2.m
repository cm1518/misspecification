%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example file for composite likelihood with two univariate models: AR(1) and MA(1)
%where the true dgp is either an AR(2) or an MA(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%options for simulating data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model=2; %=1:AR(2), =2: MA(2)

%sample size
T=200;

%joint parameter

sigma=1;

%AR(1) parameter
rho1=.5;
%AR(2) parameter
rho2=0.2;

%MA(1) parameter
theta1=.7;
%MA(2) parameter
theta2=.1;


data_temp=[0;0];
eps_lagged=[0;0];
if model==1

for t=1:T
data_vec(t)=rho1*data_temp(1)+rho2*data_temp(2)+sigma*randn;
data_temp=[data_vec(t);data_temp(1)];
end

elseif model==2

for t=1:T
eps_current=randn;
data_temp=sigma*eps_current+theta1*eps_lagged(1)+theta2*eps_lagged(2);
eps_lagged=[eps_current;eps_lagged(1)];
data_vec(t)=data_temp;
end
end
%both models share the same observables in this example

data{1}=data_vec;
data{2}=data_vec;


save('data.mat','data');


%now set up the estimation

setup_simple_example
%and run the estimation

[ draws, acc_rate, log_posteriors, statedraws, individual_post_kernels] = sampling_MH( setup );

figure;
subplot(5,1,1)
hist(draws(1,:),100)
title('\sigma')
subplot(5,1,2)
hist(draws(2,:),100)
title('\rho')
subplot(5,1,3)
hist(draws(3,:),100)
title('\theta')
subplot(5,1,4)
hist(draws(4,:),100)
title('weight on AR(1)')
subplot(5,1,5)
hist(draws(5,:),100)
title('weight on MA(1)')


