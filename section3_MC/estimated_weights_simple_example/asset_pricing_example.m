
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example file for composite likelihood with two univariate models: one for 
% dividends  and the other  for  stock  returns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the model  is
% d(t)= e(t)-theta* e(t-1) 0<theta<1
% p(t)=(1-beta *theta) e_t - theta* e(t-1)    var(e(t))=sigma
% generate beta  as  Uniform random variable with  range  [0.93, 0.97]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%options for simulating data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%model=1:dividends  2: stock  prices

%sample size
T=100;

%joint parameters
sigma=1;
theta=.5;

% stock  price parameter
bet=.95;

  uu=[0;randn(T,1)];

  betr=[.93+0.04*rand(T,1)];   
for t=1:T
data_vec1(t)=sigma*uu(t+1,1)- theta*sigma*uu(t,1);
end

for t=1:T
%data_vec2(t)=(1-bet*theta)*sigma*uu(t+1,1)- theta*sigma*uu(t,1);
data_vec2(t)=(1-betr(t,1)*theta)*sigma*uu(t+1,1)- theta*sigma*uu(t,1);
end

 
data{1}=data_vec1;
data{2}=data_vec2;


save('data.mat','data');


%now set up the estimation

setup_asset_pricing;

%and run the estimation
[ draws, acc_rate, log_posteriors, statedraws, individual_post_kernels] = sampling_MH( setup );


figure;
subplot(4,1,1)
hist(draws(1,:),100)
title('\sigma')
subplot(4,1,2)
hist(draws(2,:),100)
title('\theta')
subplot(4,1,3)
hist(draws(3,:),100)
title('\beta')
subplot(4,1,4)
hist(draws(4,:),100)
title('weight on dividend model')
