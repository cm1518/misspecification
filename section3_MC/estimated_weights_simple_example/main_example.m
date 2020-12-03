%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example file for composite likelihood with two univariate models: AR(1) and MA(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%options for simulating data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model=1; %=1:AR(1), =2: MA(1)

%sample size
T=10;

%joint parameter

sigma=1;

%AR(1) parameter
rho=.5;

%MA(1) parameter
theta=.5;

data_temp=0;
eps_lagged=0;
if model==1

for t=1:T
data_temp=rho*data_temp+sigma*randn;
data_vec(t)=data_temp;
end

elseif model==2

for t=1:T
eps_current=randn;
data_temp=sigma*eps_current+theta*eps_lagged;
eps_lagged=eps_current;
data_vec(t)=data_temp;
end
end
%both models shMAe the same observables in this example

data{1}=data_vec;
data{2}=data_vec;


save('data.mat','data');


%now set up the estimation

setup_CL

%and run the estimation

[ draws, acc_rate, log_posteriors, statedraws, individual_post_kernels] = sampling_MH( setup );

figure;
subplot(4,1,1)
hist(draws(1,:),100)
title('\sigma')
subplot(4,1,2)
hist(draws(2,:),100)
title('\rho')
subplot(4,1,3)
hist(draws(3,:),100)
title('\theta')
subplot(4,1,4)
hist(draws(4,:),100)
title('weight on AR(1) model')


figure;
plot(squeeze(individual_post_kernels)')
title('posterior kernels of each model')
legend('AR(1)','MA(1)')

%compute marginal data density for CL model

marginal_dens_full = m_harmonic(draws',log_posteriors);

