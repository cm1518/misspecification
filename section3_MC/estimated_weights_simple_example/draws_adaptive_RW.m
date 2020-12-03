function [ acc_rate ,param_draws,log_posteriors,x_draws,add_matrices_draws] = draws_adaptive_RW( posterior_draw,cov_matrix,old_posterior,setup,data ,mean )
%returns posterior draws for the standard MH algorithm with a RW propsal
%density

acceptances=0;
param_draws=zeros(setup.length_param_vector,setup.number_of_draws/setup.keep_draw);
temp_draw=inv_transform( posterior_draw,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );   


log_posteriors=zeros(1,setup.number_of_draws/setup.keep_draw);

x_draws=zeros(setup.state_size,setup.sample_size+1,setup.number_of_draws/setup.keep_draw);

if setup.add_matrices
add_matrices_draws=zeros(setup.dim_add_matrices(1),setup.dim_add_matrices(2),setup.number_of_draws/setup.keep_draw);
else
 add_matrices_draws=[];   
end


    
%    if setup.likelihood==1
%       % ll_function=SSKF_wrap( params,setup,data );
%        posterior=@(params,setup,data) (prior(params,setup)+SSKF_wrap( params,setup,data ));
%        
%    elseif setup.likelihood==2
%        %ll_function=KF_wrap;
%         posterior=@(params,setup,data) (prior(params,setup)+KF_wrap( params,setup,data ));
%      
%    else
%   
%    end
size(cov_matrix)
size(setup.scaling_adaptive*setup.eps_adaptive*eye(setup.length_param_vector))
[trash xdraw add_matrix]=posteriorwithstates(posterior_draw,setup,data);
cov_initial=setup.scaling_adaptive*cov_matrix+setup.scaling_adaptive*setup.eps_adaptive*eye(setup.length_param_vector); %initial cov matrix for adaptation
for nn=1:setup.number_of_draws
   if (nn/setup.disp_iter)==floor(nn/setup.disp_iter)
    nn
    acceptances/nn
   % diag(cov_initial)
   end
   param_prop=proposal_draw(posterior_draw,cov_initial, 1,setup); %setting scaling to 1 since it is now incorporated in cov_initial
param_prop=param_prop';
[post_prop x_prop add_matrix_prop]=posteriorwithstates(param_prop,setup,data);
%post_prop=posterior(param_prop,setup,data);
alpha=min(1,exp(post_prop-old_posterior));
if rand<alpha
   posterior_draw=param_prop;
   old_posterior=post_prop;
   acceptances=acceptances+1;
   temp_draw=inv_transform( posterior_draw,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );
   xdraw=x_prop;
   add_matrix=add_matrix_prop;
end

[cov_initial,mean]=recursive_moments(nn+setup.burn_in, posterior_draw,cov_initial, mean,setup );


if ((nn)/setup.keep_draw)==floor(((nn)/setup.keep_draw))
param_draws(:,(nn)/setup.keep_draw)=temp_draw;
log_posteriors((nn)/setup.keep_draw)=old_posterior;
x_draws(:,:,(nn)/setup.keep_draw)=xdraw;


if setup.add_matrices
add_matrices_draws(:,:,(nn)/setup.keep_draw)=add_matrix;
end

end

end


acc_rate=acceptances/setup.number_of_draws;

end
