function [ log_l, indicator, epsilon, uvec] = likelihood(data,  params,setup)


%unwrapping parameters
[ params ] = params_mod( params,setup );
mean_ind=mean(setup.indicator);

temp_sum1=sum(params.alpha_gen_neg.*repmat((1+(params.beta_gen_neg*(min(setup.indicator)-mean_ind))),1,size(params.alpha_gen_neg,2)),2);
temp_sum2=sum(params.alpha_gen_neg.*repmat((1+(params.beta_gen_neg*(max(setup.indicator)-mean_ind))),1,size(params.alpha_gen_neg,2)),2);

temp_sum3=sum(params.alpha_gen_pos.*repmat((1+(params.beta_gen_pos*(min(setup.indicator)-mean_ind))),1,size(params.alpha_gen_pos,2)),2);
temp_sum4=sum(params.alpha_gen_pos.*repmat((1+(params.beta_gen_pos*(max(setup.indicator)-mean_ind))),1,size(params.alpha_gen_pos,2)),2);


if min([temp_sum1(setup.index_unrestricted) temp_sum2(setup.index_unrestricted) temp_sum3(setup.index_unrestricted) temp_sum4(setup.index_unrestricted)])<0
log_l=-inf;
indicator=[];
epsilon=[];
uvec=[];
else

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

%options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',500, 'MaxFunEvals',5000);

%build X matrix of RHS trend variables
X=[];
for jj=1:setup.polynomials+1
   % X=[X [0:setup.sample_size-1]'.^(jj-1)];
    X=[X (2*[0:setup.sample_size-1]'/(setup.sample_size-1)-1).^(jj-1)]; %rescale polynomials to make algorithm more stable
    % see https://www.rssd.esa.int/SP/LISAPATHFINDER/docs/Data_Analysis/DA_Nine/Annex7b.pdf
end




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
[ Sigma] = unwrap_NL_lagged2(  params.alpha_gen_neg,params.beta_gen_neg,params.b_gen_neg,params.c_gen_neg,params.alpha_gen_pos,params.beta_gen_pos,params.b_gen_pos,params.c_gen_pos,fliplr(eps_initial),setup,setup.indicator(tt:tt+setup.lags-1),mean_ind);






condmean=zeros(setup.size_obs,1);

for ll=1:lag_length
   
   
   condmean=condmean+Sigma(:,:,ll)*eps_initial(:,end-ll); %note the reverse ordering of Sigma and epsilon
end
condmean=condmean+[params.intercepts params.poly_coefficients]*X(tt,:)';
u=data(:,tt)-condmean;
uvec(:,tt)=u;
temp_neg=setup.store_responses(:,:,1);
temp2=sum(params.alpha_gen_neg.*repmat((1+(params.beta_gen_neg*(setup.indicator(tt+setup.lags)-mean_ind))),1,size(params.alpha_gen_neg,2)),2);
temp_neg(setup.index_unrestricted:end,setup.index_unrestricted)=params.r_neg.*[temp2(setup.index_unrestricted:end,:)];



epsilon_temp=temp_neg\u;

if epsilon_temp(setup.index_unrestricted)>0
     ind=1;
  temp_pos=setup.store_responses(:,:,1);
 temp2=sum(params.alpha_gen_pos.*repmat((1+(params.beta_gen_pos*(setup.indicator(tt+setup.lags)-mean_ind))),1,size(params.alpha_gen_pos,2)),2);

temp_pos(setup.index_unrestricted:end,setup.index_unrestricted)=params.r_pos.*[temp2(setup.index_unrestricted:end,:)];
temp_contemp=temp_pos;
epsilon_temp=temp_pos\u;
else
    ind=0;
   temp_contemp=temp_neg;
end


Sigma_contemp=temp_contemp;
covariance=Sigma_contemp*Sigma_contemp';
epsilon(:,tt+lag_length)=epsilon_temp;
%log_l_temp_vec=zeros(setup.size_obs,1);
%for jj=1:setup.size_obs
log_l_temp=-(setup.size_obs/2)*log(2*pi)-.5*log(det(covariance))-.5*((Sigma_contemp*epsilon_temp)'/covariance)*(Sigma_contemp*epsilon_temp);
%end



log_l=log_l+log_l_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
if (sum(epsilon(3,:)==0)-setup.lags)>2
    log_l=-inf;
end
end