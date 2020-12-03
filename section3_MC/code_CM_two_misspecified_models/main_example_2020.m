% main_example_2020  final revision 13-10-2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% composite likelihood estimation with univariate AR(1) and MA(1) models
% true dgp is either an ARMA(1), an AR(2), or an MA(2)
% compute  cl with random omega, equal omega, and  MSE omega compare with  
% with ar(1), ma(1) and ARMA(1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collects MSE and  KL  for  various  models
% computes posterior  omega  and  BMA weight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: prior  for  estimation needs  to  be  changed  depending  on  the 
% parameters of  the  DGP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc;
clear;
rng('default')
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options for simulating data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model=3; %=1:AR(2), =2: MA(2), =3: ARMA(1,1)
T=250;  % sampe size

% common parameter
sigma=1.0; %0.2 %0.5; %0.8 %1.0; %1.2 

%AR(1) parameter
rho1=0.3; %0.3; %0.9; %0.6;
%AR(2) parameter
rho2=0.0;

%MA(1) parameter
theta1=0.8; %0.8; %0.2; %0.5;
%MA(2) parameter
theta2=0.0;

data_temp=[0;0];
eps_lagged=[0;0];

% random number  generator
abc=randn(T,1);
%data_vec=zeros(1,T);
if model==1
    for t=1:T
      data_vec(t)=rho1*data_temp(1)+rho2*data_temp(2)+sigma*abc(t,1);
      data_temp=[data_vec(t);data_temp(1)];
    end

elseif model==2
   for t=1:T
     eps_current=abc(t,1); %randn;
     data_temp=sigma*eps_current+theta1*eps_lagged(1)+theta2*eps_lagged(2);
     eps_lagged=[eps_current;eps_lagged(1)];
     data_vec(t)=data_temp;
   end

elseif model==3
   for t=1:T
     eps_current=abc(t,1); %randn;
     data_vec(t)=rho1*data_temp(1)+ sigma*eps_current+theta1*eps_lagged(1);
     eps_lagged=[eps_current;eps_lagged(1)];
     data_temp=[data_vec(t);data_temp(1)];
   end
end


% both models share the same observables
data{1}=(data_vec);
data{2}=(data_vec);

save('data.mat','data');

cd ../estimated_weights_simple_example
save('data.mat','data');
data=data_vec;
cd ../single_model
save data data

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% 1) CL with estimated weights     %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

cd ../estimated_weights_simple_example

% parameters of  the  MCMC exercise
clear setup
setup_simple_example;
setup.sample_size=T;
setup.number_of_draws=50000;
setup.scaling_draws=5000;
setup.check_scaling=200;
setup.keep_draw=10;


% for ARMA specification with AR=0.6 and MA=0.5 use
% setup.normal_prior_means=[0;0.6;0]; %for prior means
% setup.normal_prior_std=[.18;.2;.2]; % prior std [0.15, .2,.2];

% for ARMA specification with AR=0.6 and MA=0.8 use
%setup.normal_prior_means=[0;0.5;0.0]; %for prior means
%setup.normal_prior_std=[.15; .2; .2];

% for ARMA specification with AR=0.9 and MA=0.2 use
% setup.normal_prior_means=[0;0.6;0]; %for prior means
% setup.normal_prior_std=[.18;.2;.2]; % prior std [0.15, .2,.2];

% for ARMA specification  with AR=0.3and M=0.8 use
 setup.normal_prior_means=[0;0.0;0.6];  % [0, 0, 0.5] for prior means
 setup.normal_prior_std=[.18;.2;.2]; % prior std [0.15, .2,.2];

[draws_ew, acc_rate_ew, log_posteriors_ew, statedraws_ew, individual_post_kernels_ew] = sampling_MH( setup );

cd ../code_CM_two_misspecified_models

% collecting  estimates  of  distribution of  parameters
sig_ew(1)=mean(draws_ew(1,:));
sig_ew(2)=std(draws_ew(1,:));
sig_ew(3)=median(draws_ew(1,:));
sig_ew(4)=mode(draws_ew(1,:));

rho_ew(1)=mean(draws_ew(2,:));
rho_ew(2)=std(draws_ew(2,:));
rho_ew(3)=median(draws_ew(2,:));
rho_ew(4)=mode(draws_ew(2,:));

thet_ew(1)=mean(draws_ew(3,:));
thet_ew(2)=std(draws_ew(3,:));
thet_ew(3)=median(draws_ew(3,:));
thet_ew(4)=mode(draws_ew(3,:));

ome_ew(1)=mean(draws_ew(4,:));
ome_ew(2)=std(draws_ew(4,:));
ome_ew(3)=median(draws_ew(4,:));
ome_ew(4)=mode(draws_ew(4,:));

% prior draws for  omega
prior_draw=zeros(length(log_posteriors_ew),2);
for jj=1:length(log_posteriors_ew)
   prior_draw(jj,:)= drchrnd(setup.dirichlet_prior_parameters,1);     
end

% prior  and  posterior  weights
%{
figure(101);
histogram(prior_draw(:,1),100); % 'Normalization', 'pdf'); 
hold on;
histogram(draws_ew(4,:),100); %'Normalization', 'pdf'); 
hold off;
title('Weight on AR(1) model, true model ARMA(1,1)')
grid on
%}

% MSE of  the  parameters
MSE_sig=mean((draws_ew(1,:)-sigma).^2);
MSE_rho=mean((draws_ew(2,:)-rho1).^2);
MSE_thet=mean((draws_ew(3,:)-theta1).^2);

% marginal  likelihood
marginal_dens_full_ew = m_harmonic(draws_ew',log_posteriors_ew);


% KL integrating  over  parameters

% edges  of  the  histogram  bins
zx=prctile(exp(data),[10 20 30 40 50 60 70 80 90]);
simar=zeros(length(draws_ew),T); simma=zeros(length(draws_ew),T);
simcl=zeros(length(draws_ew),T);
zkl=zeros(length(draws_ew),3,4);
d1=0.0; d2=0.0; d3=0.0;
for  kk=1:length(draws_ew)

    % simulating  from  an  AR(1)
data_temp=[0;0]; %dataAR=zeros(T,1);  
for t=1:T
dataAR(t)=draws_ew(2,kk)*data_temp(1)+draws_ew(1,kk)*abc(t,1);
data_temp=[dataAR(t);data_temp(1)];
end
simar(kk,:)=exp(dataAR(1,:));

% simulating data  from  the  MA1
eps_lagged=[0;0]; %dataMA=zeros(T,1);
for t=1:T
eps_current=abc(t,1); %randn;
data_temp=draws_ew(1,kk)*eps_current+draws_ew(3,kk)*eps_lagged(1);
eps_lagged=[eps_current;eps_lagged(1)];
dataMA(t)=data_temp;
end
simma(kk,:)=exp(dataMA(1,:));

% simulating from  CL
%simcl(kk,:)=draws_ew(4,kk)*simar(kk,:)+(1-draws_ew(4,kk))*simma(kk,:);  
simcl(kk,:)=((simar(kk,:)).^(draws_ew(4,kk))).*((simma(kk,:)).^(1-draws_ew(4,kk)));

aa=hist(exp(data),zx);      % histogram of simulated  data
ww=hist(simar(kk,:),zx);    % histogram of  ar(1)
vv=hist(simma(kk,:),zx);    % histogram of  ma(1)
zz=hist(simcl(kk,:),zx);    % histogram of cl

% KL divergence % 
akl3=0.0; akl1=0.0; akl2=0.0; 
for i=1:length(aa')
    if ww(i) ~=0
        if  aa(i) ~=0
        akl1=akl1-aa(i)*log(ww(i)/aa(i)); % ar(1) divergence
        end
    end
    if  vv(i)~=0
        if aa(i) ~=0
        akl2=akl2-aa(i)*log(vv(i)/aa(i)); % ma(1) divergence
        end
    end
    if zz(i) ~=0
        if aa(i) ~=0
        akl3=akl3-aa(i)*log(zz(i)/aa(i)); % CL divergence
        end
    end
    end
zkl(kk,1,1)=akl1;zkl(kk,2,1)=akl2; zkl(kk,3,1)=akl3;

d1=d1+(1.0/length(draws_ew))*akl1;
d2=d2+(1.0/length(draws_ew))*akl2;
d3=d3+(1.0/length(draws_ew))*akl3;

end

 


% return


 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% 2) CL estimation with equal weights  %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

cd ../code_CM_two_misspecified_models

% parameters of  the  MCMC exercise
setup_CL
setup.number_of_draws=50000;
setup.scaling_draws=5000;
setup.check_scaling=200;
setup.keep_draw=10;


%For  ar=0.6, ma=0.5
%setup.normal_prior_means=[0;0.6;0]; %prior means
%setup.normal_prior_std=[.18;.2;.2]; % prior std [.15,.2,.2]

%For arma specification with  ar=0.6, ma=0.8
%setup.normal_prior_means=[0;0.5;0]; %prior means
%setup.normal_prior_std=[.15;.2;.2]; % prior std [.15,.2,.2]


% for ARMA specification with AR=0.9 and MA=0.2 use
% setup.normal_prior_means=[0;0.6;0]; %for prior means
%setup.normal_prior_std=[.18;.2;.2]; % prior std [0.15, .2,.2];

% for ARMA specification  with AR=0.3and M=0.8 use
 setup.normal_prior_means=[0;0.0;0.6];  %for prior means
 setup.normal_prior_std=[.18;.2;.2]; % prior std [0.15, .2,.2];

[ draws, acc_rate, log_posteriors, statedraws, individual_post_kernels] = sampling_MH( setup );

% statistics  of  estimated  distribution
sig_55(1)=mean(draws(1,:));
sig_55(2)=std(draws(1,:));
sig_55(3)=median(draws(1,:));
sig_55(4)=mode(draws(1,:));

rho_55(1)=mean(draws(2,:));
rho_55(2)=std(draws(2,:));
rho_55(3)=median(draws(2,:));
rho_55(4)=mode(draws(2,:));

thet_55(1)=mean(draws(3,:));
thet_55(2)=std(draws(3,:));
thet_55(3)=median(draws(3,:));
thet_55(4)=mode(draws(3,:));

 
 
% marginal  likelihood
marginal_dens_55 = m_harmonic(draws',log_posteriors);

% MSE of  the  parameters
MSE_sig55=mean((draws(1,:)-sigma).^2);
MSE_rho55=mean((draws(2,:)-rho1).^2);
MSE_thet55=mean((draws(3,:)-theta1).^2);


  
 
% kl integrating  over  parameter draws
d155=0.0; d255=0.0; d355=0.0;
simar55=zeros(length(draws_ew),T); simma55=zeros(length(draws_ew),T);
simcl55=zeros(length(draws_ew),T);
for  kk=1:length(draws)
%dataAR55=zeros(T,1); dataMA55=zeros(T,1);
data_temp=[0;0];
    for t=1:T
        dataAR55(t)=draws(2,kk)*data_temp(1)+draws(1,kk)*abc(t,1);
        data_temp=[dataAR55(t);data_temp(1)];
    end
simar55(kk,:)=exp(dataAR55(1,:));

eps_lagged=[0;0];
    for t=1:T
        eps_current=abc(t,1); %randn;
        data_temp=draws(1,kk)*eps_current+draws(3,kk)*eps_lagged(1);
        eps_lagged=[eps_current;eps_lagged(1)];
        dataMA55(t)=data_temp;
    end
simma55(kk,:)=exp(dataMA55(1,:));

%simcl55(kk ,:)=0.5*simar55(kk,:)+0.5*simma55(kk,:);    
simcl55(kk,:)=((simar55(kk,:)).^(0.5)).*((simma55(kk,:)).^(0.5));

ww=hist(simar55(kk,:),zx);    % histogram of  simulated ar(1)
vv=hist(simma55(kk,:),zx);    % histogram of  simulated ma(1)
zz=hist(simcl55(kk,:),zx);      % histogram of cl

% KL divergence % 
akl3=0.0; akl1=0.0; akl2=0.0; 
for i=1:length(aa')
    if ww(i) ~=0
        if  aa(i) ~=0
        akl1=akl1-aa(i)*log(ww(i)/aa(i)); % ar(1) divergence
        end
    end
    if  vv(i)~=0
        if aa(i) ~=0
        akl2=akl2-aa(i)*log(vv(i)/aa(i)); % ma(1)   divergence
        end
    end
    if zz(i) ~=0
        if aa(i) ~=0
        akl3=akl3-aa(i)*log(zz(i)/aa(i)); % cl divergence
        end
    end
end
zkl(kk,1,2)=akl1; zkl(kk,2,2)=akl2;  zkl(kk,3,2)=akl3;

d155=d155+(1.0/length(draws))*akl1;
d255=d255+(1.0/length(draws))*akl2;
d355=d355+(1.0/length(draws))*akl3;

end

 
 
 
 
 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% 3) CL estimation with  mse weights %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% with AR(2)
%setup.model_weights=[.98 .02]; %model weights

% with ARMA
setup.number_of_draws=50000;
setup.scaling_draws=5000;
setup.check_scaling=200;
setup.keep_draw=10;

% AR(0.6), MA(0.5)
%setup.model_weights=[.97 .03]; 
%setup.normal_prior_means=[0;0.6;0]; %prior means
%setup.normal_prior_std=[.18;.2;.2]; % prior std  [.15,.2.,2]


% AR(0.6), MA(0.8)
%setup.model_weights=[.97 .03]; 
%setup.normal_prior_means=[0;0.5;0]; %prior means
%setup.normal_prior_std=[.2;.2;.2]; % prior std  [.15,.2.,2]

% AR(0.9), MA(0.2)
%setup.model_weights=[.98 .02]; 
%setup.normal_prior_means=[0;0.6;0]; %prior means
%setup.normal_prior_std=[.18;.2;.2]; % prior std  [.15,.2.,2]

% AR(0.3), MA(0.8)
setup.model_weights=[.94 .06]; 
setup.normal_prior_means=[0;0.0;0.6]; %prior means
setup.normal_prior_std=[.18;.2;.2]; % prior std  [.15,.2.,2]



%setup.normal_prior_means=[0;0.6;0.6]; %prior means


[ drawsmse, acc_ratemse, log_posteriorsmse, statedrawsmse, individual_post_kernelsmse] = sampling_MH( setup );

% statistics  of  posterior
sig_mse(1)=mean(drawsmse(1,:));
sig_mse(2)=std(drawsmse(1,:));
sig_mse(3)=median(drawsmse(1,:));
sig_mse(4)=mode(drawsmse(1,:));

rho_mse(1)=mean(drawsmse(2,:));
rho_mse(2)=std(drawsmse(2,:));
rho_mse(3)=median(drawsmse(2,:));
rho_mse(4)=mode(drawsmse(2,:));

thet_mse(1)=mean(drawsmse(3,:));
thet_mse(2)=std(drawsmse(3,:));
thet_mse(3)=median(drawsmse(3,:));
thet_mse(4)=mode(drawsmse(3,:));

 
 
% marginal  likelihood
marginal_dens_mse = m_harmonic(drawsmse',log_posteriorsmse);

% MSE of  the  parameters
MSE_sigmse=mean((drawsmse(1,:)-sigma).^2);
MSE_rhomse=mean((drawsmse(2,:)-rho1).^2);
MSE_thetmse=mean((drawsmse(3,:)-theta1).^2);


 
% KL integrating  over  parameters
d1mse=0.0; d2mse=0.0; d3mse=0.0;
simarmse=zeros(length(draws_ew),T); simmamse=zeros(length(draws_ew),T);
simclmse=zeros(length(draws_ew),T);

for  kk=1:length(drawsmse)
% simuating form  an  AR(1)
    data_temp=[0;0]; %dataARmse=zeros(T,1); dataMAmse=zeros(T,1);
for t=1:T
dataARmse(t)=drawsmse(2,kk)*data_temp(1)+drawsmse(1,kk)*abc(t,1);
data_temp=[dataARmse(t);data_temp(1)];
end
simarmse(kk,:)=exp(dataARmse(1,:));

% simulating data  from  the  MA1
eps_lagged=[0;0];
for t=1:T
eps_current=abc(t,1); %randn;
data_temp=drawsmse(1,kk)*eps_current+drawsmse(3,kk)*eps_lagged(1);
eps_lagged=[eps_current;eps_lagged(1)];
dataMAmse(t)=data_temp;
end
simmamse(kk,:)=exp(dataMAmse(1,:));
% simulating form  cl
%simclmse(kk,:)=setup.model_weights(1)*simarmse(kk,:)+setup.model_weights(2)*simmamse(kk,:);  %linear
simclmse(kk,:)=((simarmse(kk,:)).^(setup.model_weights(1))).*((simmamse(kk,:)).^(setup.model_weights(2)));

ww=hist(simarmse(kk,:),zx);    % histogram of  simulated ar(1)
vv=hist(simmamse(kk,:),zx);    % histogram of  simulated ma(1)
zz=hist(simclmse(kk,:),zx);      % histogram of cl

% KL divergence % 
akl3=0.0; akl1=0.0; akl2=0.0; 
for i=1:length(aa')
    if ww(i) ~=0
        if  aa(i) ~=0
        akl1=akl1-aa(i)*log(ww(i)/aa(i)); % ar(1) divergence
        end
    end
    if  vv(i)~=0
        if aa(i) ~=0
        akl2=akl2-aa(i)*log(vv(i)/aa(i)); % ma(1)   divergence
        end
    end
    if zz(i) ~=0
        if aa(i) ~=0
        akl3=akl3-aa(i)*log(zz(i)/aa(i)); % cl divergence
        end
    end
end
zkl(kk,1,3)=akl1; zkl(kk,2,3)=akl2; zkl(kk,3,3)=akl3;

d1mse=d1mse+(1.0/length(drawsmse))*akl1;
d2mse=d2mse+(1.0/length(drawsmse))*akl2;
d3mse=d3mse+(1.0/length(drawsmse))*akl3;
end

 
   

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% 4) estimation of AR(1) %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')

cd ../single_model

% options  for  MCMC
setup_AR;
setup.sample_size=T;
setup.number_of_draws=50000;
setup.scaling_draws=5000;
setup.check_scaling=200;
setup.keep_draw=10;
setup.normal_prior_means=[0;.6]; %prior means
setup.normal_prior_std=[.2;.2]; % prior std
[ drawsAR, acc_rateAR, log_posteriorsAR, statedrawsAR, individual_post_kernelsAR] = sampling_MH( setup );

cd ../code_CM_two_misspecified_models

% statistics of posterior
sig_ar(1)=mean(drawsAR(1,:));
sig_ar(2)=std(drawsAR(1,:));
sig_ar(3)=median(drawsAR(1,:));
sig_ar(4)=mode(drawsAR(1,:));

rho_ar(1)=mean(drawsAR(2,:));
rho_ar(2)=std(drawsAR(2,:));
rho_ar(3)=median(drawsAR(2,:));
rho_ar(4)=mode(drawsAR(2,:));

% marginal  likelihood 
marginal_dens_AR = m_harmonic(drawsAR',log_posteriorsAR);

% MSE of  the  parameters
MSE_sigAR=mean((drawsAR(1,:)-sigma).^2);
MSE_rhoAR=mean((drawsAR(2,:)-rho1).^2);

 
% KL integrating  over  parameters
d1AR=0.0;
simarAR=zeros(length(draws_ew),T); 
for  kk=1:length(drawsAR)
    data_temp=[0;0];%dataARAR=zeros(T,1); 
        for t=1:T
            dataARAR(t)=drawsAR(2,kk)*data_temp(1)+drawsAR(1,kk)*abc(t,1);
            data_temp=[dataARAR(t);data_temp(1)];
        end
    simarAR(kk,:)=exp(dataARAR(1,:));
    
ww=hist(simarAR(kk,:),zx);    % histogram of  simulated ar(1)
 
% KL divergence % 
akl1=0.0; 
for i=1:length(aa')
     if ww(i) ~=0
        if  aa(i) ~=0
        akl1=akl1-aa(i)*log(ww(i)/aa(i)); % ma(1) divergence
        end
     end
end
zkl(kk,1,4)=akl1;
d1AR=d1AR+(1.0/length(drawsAR))*akl1;
end


 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% 5) MA(1) estimation    %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')

cd ../single_model

% parameters of  MCMC
setup_MA;
setup.sample_size=T;
setup.number_of_draws=50000;
setup.scaling_draws=5000;
setup.check_scaling=200;
setup.keep_draw=10;
setup.normal_prior_means=[0;0]; %prior means
setup.normal_prior_std=[.2;.2]; % prior std
[ drawsMA, acc_rateMA, log_posteriorsMA, statedrawsMA, individual_post_kernelsMA] = sampling_MH( setup );

cd ../code_CM_two_misspecified_models

% statistics  of  posterior
sig_ma(1)=mean(drawsMA(1,:));
sig_ma(2)=std(drawsMA(1,:));
sig_ma(3)=median(drawsMA(1,:));
sig_ma(4)=mode(drawsMA(1,:));

thet_ma(1)=mean(drawsMA(2,:));
thet_ma(2)=std(drawsMA(2,:));
thet_ma(3)=median(drawsMA(2,:));
thet_ma(4)=mode(drawsMA(2,:));

% marginal  likelihood 
marginal_dens_MA = m_harmonic(drawsMA',log_posteriorsMA);

% MSE of  the  parameters
MSE_sigMA=mean((drawsMA(1,:)-sigma).^2);
MSE_thetMA=mean((drawsMA(2,:)-theta1).^2);

 
 % KL integrating  over  parameters
 d1MA=0.0; simmaMA=zeros(length(draws_ew),T); 
for  kk=1:length(drawsMA)
eps_lagged=[0;0]; %dataMAMA=zeros(T,1);
    for t=1:T
        eps_current=abc(t,1); %randn;
        data_temp=drawsMA(1,kk)*eps_current+drawsMA(2,kk)*eps_lagged(1);
        eps_lagged=[eps_current;eps_lagged(1)];
        dataMAMA(t)=data_temp;
    end
simmaMA(kk,:)=exp(dataMAMA(1,:));
vv=hist(simmaMA(kk,:),zx);    % histogram of  simulated ma(1)

 
 akl1=0.0; 
for i=1:length(aa')
     if  vv(i) ~=0
        if aa(i) ~=0
            akl1=akl1-aa(i)*log(vv(i)/aa(i)); % cl   divergence
        end
     end
end
zkl(kk,2,4)=akl1;
d1MA=d1MA+(1.0/length(drawsMA))*akl1;
end



 
 

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% 6) ARMA(1,1) estimation   %%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

cd ../single_model
% parametes of  MCMC

setup_ARMA;
setup.sample_size=T;
setup.number_of_draws=50000;
setup.scaling_draws=5000;
setup.check_scaling=200;
setup.keep_draw=10;
setup.initial_parameter=[1;.5;.1]; % [.1;.8;.1]
setup.normal_prior_means=[0;.6;0]; %prior means
setup.normal_prior_std=[.2;.2;.2]; % prior std
[ drawsARMA, acc_rateARMA, log_posteriorsARMA, statedrawsARMA, individual_post_kernelsARMA] = sampling_MH( setup );

cd ../code_CM_two_misspecified_models
% statistics  of  posterior
sig_arma(1)=mean(drawsARMA(1,:));
sig_arma(2)=std(drawsARMA(1,:));
sig_arma(3)=median(drawsARMA(1,:));
sig_arma(4)=mode(drawsARMA(1,:));

rho_arma(1)=mean(drawsARMA(2,:));
rho_arma(2)=std(drawsARMA(2,:));
rho_arma(3)=median(drawsARMA(2,:));
rho_arma(4)=mode(drawsARMA(2,:));

thet_arma(1)=mean(drawsARMA(3,:));
thet_arma(2)=std(drawsARMA(3,:));
thet_arma(3)=median(drawsARMA(3,:));
thet_arma(4)=mode(drawsARMA(3,:));

%compute likleihood
marginal_dens_ARMA = m_harmonic(drawsARMA',log_posteriorsARMA);

% MSE of  the  parameters
MSE_sigARMA=mean((drawsARMA(1,:)-sigma).^2);
MSE_rhoARMA=mean((drawsARMA(2,:)-rho1).^2);
MSE_thetARMA=mean((drawsARMA(3,:)-theta1).^2);


  
% KL integrating  over  parameters
d1ARMA=0.0;simARMA=zeros(length(draws_ew),T);
for  kk=1:length(drawsARMA)
    data_temp=[0;0];  eps_lagged=[0;0];
        for t=1:T
        eps_current=abc(t,1); %randn;
        data_vec(1,t)=drawsARMA(2,kk)*data_temp(1)+ ...
             drawsARMA(1,kk)*eps_current+drawsARMA(3,kk)*eps_lagged(1);
        eps_lagged=[eps_current;eps_lagged(1)];
        data_temp=[data_vec(1,t);data_temp(1)];
        end
 simARMA(kk,:)=exp(data_vec(1,:));

ww=hist(simARMA(kk,:),zx);    % histogram of  simulated arma(1,1)
 
% KL divergence % 
akl1=0.0; 
for i=1:length(aa')
      if ww(i) ~=0
        if  aa(i) ~=0
        akl1=akl1-aa(i)*log(ww(i)/aa(i)); % ma(1) divergence
        end
      end
end
zkl(kk,3,4)=akl1;
d1ARMA=d1ARMA+(1.0/length(drawsARMA))*akl1;
end
 
 



 % eliminating  average kl  measures  which are  negative
zklm=zeros(1,size(zkl,2),4); num=zeros(size(zkl,2),4);
akl=zeros(2,size(zkl,2),4);
for m=1:4
 for jj=1:size(zkl,2)
     for i=1:length(draws)
       if zkl(i,jj,m)>0
         zklm(1,jj,m)=zklm(1,jj,m)+zkl(i,jj,m);
         num(jj,m)=num(jj,m)+1;
       end
     end
   akl(1,jj,m)=zklm(1,jj,m)/num(jj,m);
   akl(2,jj,m)=median(zkl(:,jj,m));  % taking the  median  value across draws
end
end


%clc;

 
disp('MSE of  sigma')
disp([MSE_sig MSE_sig55 MSE_sigmse MSE_sigAR MSE_sigMA])

disp('average predictive  KL measures')
rkl(1,1)=akl(1,3,1);
rkl(1,2)=akl(2,3,1);
rkl(2,1)=akl(1,3,2);
rkl(2,2)=akl(2,3,2);
rkl(3,1)=akl(1,3,3);
rkl(3,2)=akl(2,3,3);
rkl(4,1)=akl(1,1,4);
rkl(4,2)=akl(2,1,4);
rkl(5,1)=akl(1,2,4);
rkl(5,2)=akl(2,2,4);

disp(rkl)

 
disp('BMA weight  on  AR(1)')
disp(exp(marginal_dens_AR)/(exp(marginal_dens_AR)+exp(marginal_dens_MA)))

disp('omega weight  on AR(1) for CLEsW')
disp([mode(draws_ew(4,:)), std(draws_ew(4,:))])

 
toc;
 
%save ARMA_100308_T250  