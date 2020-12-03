function [ log_l, indicator, epsilon, uvec] = likelihood(data,  params,setup)
%computes the log likelihood for the MA model
 eps_initial =setup.initial_eps; 
lag_length=setup.lags;


indicator=zeros(size(data,1),lag_length+size(data,2));
epsilon=zeros(size(data,1),lag_length+size(data,2));
%data is supposed to be a matrix with each column containing the
%observations at one point in time
uvec=zeros(size(data,1),size(data,2));

%careful about direction of epsilon/data!
%pass params as an argument to this function
%update epsilon! DONE

%setting options for optimization

options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',500, 'MaxFunEvals',5000);


epsilon(:,1:lag_length)=eps_initial;
%indicator equal to 0 if epsilon is negative or zero, 1 otherwise
log_l=0;
for tt=1:size(data,2)
%     if tt==1
%         epsilon(:,tt:tt+lag_length-1)
%     end

%solve for contemporaneous epsilons

eps_initial=[epsilon(:,tt:tt+lag_length-1) zeros(setup.size_obs,1)]; % we pad zeros on for the contmeporaneous epsilons

%first, get lagged sigma matrices
[ Sigma, intercept,alpha_diag,beta_diag] = unwrap_NL_lagged( params,fliplr(eps_initial),setup );%note the reverse ordering of Sigma and epsilon






condmean=zeros(setup.size_obs,1);

for ll=1:lag_length
   
   
   condmean=condmean+Sigma(:,:,ll)*eps_initial(:,end-ll+1); %note the reverse ordering of Sigma and epsilon
end
condmean=condmean+intercept;
u=data(:,tt)-condmean;
uvec(:,tt)=u;
epsilon_temp=zeros(size(data,1),1);

for kk=1:setup.size_obs
    [ sum_temp] = within_period_NL( params,epsilon_temp,setup,kk );
    
postopt2=@(shocks) NLcontemp_abs(shocks,alpha_diag(kk),beta_diag(kk),condmean(kk)-sum_temp);
  

[epsilon_temp(kk),~]=fzero(postopt2,0,options);



end

% [ sigma_temp] = unwrap_NL_contemp( params,epsilon_vec,setup );
% 
% covariance=sigma_temp*sigma_temp';

epsilon(:,tt+lag_length)=epsilon_temp;
log_l_temp_vec=zeros(setup.size_obs,1);
for jj=1:setup.size_obs
log_l_temp_vec(jj)=-(1/2)*log(2*pi)-.5*((epsilon_temp(jj))')*(epsilon_temp(jj))-log(2*alpha_diag(jj)*epsilon_temp(jj)^2+1)-(beta_diag(jj)+alpha_diag(jj)*epsilon_temp(jj)^2);
end



log_l=log_l+sum(log_l_temp_vec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
