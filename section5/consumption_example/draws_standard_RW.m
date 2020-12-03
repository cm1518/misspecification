function [ acc_rate ,param_draws,log_posteriors,x_draws,add_matrices_draws] = draws_standard_RW( posterior_draw,cov_matrix,scaling,old_posterior,setup,data  )
%returns posterior draws for the standard MH algorithm with a RW propsal
%density

acceptances=zeros(setup.number_blocks,1);
param_draws=zeros(setup.length_param_vector,setup.number_of_draws/setup.keep_draw);
temp_draw=inv_transform( posterior_draw,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );   


if setup.add_matrices
add_matrices_draws=zeros(setup.dim_add_matrices(1),setup.dim_add_matrices(2),setup.number_of_draws/setup.keep_draw);
else
 add_matrices_draws=[];   
end




log_posteriors=zeros(1,setup.number_of_draws/setup.keep_draw);

x_draws=zeros(setup.state_size,setup.sample_size+1,setup.number_of_draws/setup.keep_draw);

    
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

[~, xdraw, add_matrix,log_post]=posteriorwithstates(posterior_draw,setup,data);

for nn=1:setup.number_of_draws
   if (nn/setup.disp_iter)==floor(nn/setup.disp_iter)
  disp('number of draws:')
       nn
       disp('acceptance prob:')
acceptances/nn
   end
   
   for bb=1:setup.number_blocks
  

  if length(setup.index_block{bb})==setup.length_param_vector %there is only one block
  param_prop=proposal_draw(posterior_draw(1:(min(setup.weight_index)-1)),cov_matrix{bb}, scaling(bb),setup); %random walk proposal
  param_prop=param_prop';
  param_prop_dirichlet=drchrnd(setup.dirichlet_scaling*posterior_draw((min(setup.weight_index):end))/sum(posterior_draw((min(setup.weight_index):end))),1);
  %param_prop_dirichlet=drchrnd(setup.dirichlet_proposal_parameters,1);
  param_prop=[param_prop;param_prop_dirichlet];
  elseif bb<setup.number_blocks %a non-dirichlet block
  param_prop=proposal_draw(posterior_draw(setup.index_block{bb}),cov_matrix{bb}, scaling(bb),setup);
  param_prop=param_prop';
  elseif bb==setup.number_blocks %the weight block
  param_prop_dirichlet=drchrnd(setup.dirichlet_scaling*posterior_draw((min(setup.weight_index):end))/sum(posterior_draw((min(setup.weight_index):end))),1);
  param_prop=param_prop_dirichlet;
  end




 


posterior_draw_prop=posterior_draw;
posterior_draw_prop(setup.index_block{bb})=param_prop;
if length(setup.index_block{bb})==setup.length_param_vector %one block
[post_prop, x_prop, add_matrix_prop,log_post_prop]=posteriorwithstates(posterior_draw_prop,setup,data);
elseif (length(setup.index_block{bb})<setup.length_param_vector) && bb==1 %common parameter block
[post_prop, x_prop, add_matrix_prop,log_post_prop]=posteriorwithstates(posterior_draw_prop,setup,data);
elseif (length(setup.index_block{bb})<setup.length_param_vector) && bb>1 && bb<setup.number_blocks %single model block
[post_prop_mm, x_prop, add_matrix_prop2]=posteriorwithstates_single_model(posterior_draw_prop,setup,data,bb-1);
%post_prop=posterior(param_prop,setup,data);
post_prop=old_posterior-log_post(bb-1)+post_prop_mm;
elseif (length(setup.index_block{bb})<setup.length_param_vector) &&  bb==setup.number_blocks %dirichlet block
[post_prop, x_prop, add_matrix_prop3,log_post_prop]=posteriorwithstates(posterior_draw_prop,setup,data);
end




if length(setup.index_block{bb})==setup.length_param_vector
 prop_adjust=exp(dirpdf_log(posterior_draw(setup.weight_index),setup.dirichlet_scaling*posterior_draw_prop(min(setup.weight_index):end)/sum(posterior_draw_prop(min(setup.weight_index):end))))/exp(dirpdf_log(posterior_draw_prop(setup.weight_index),setup.dirichlet_scaling*posterior_draw(min(setup.weight_index):end)/sum(posterior_draw(min(setup.weight_index):end))));
   if isnan(prop_adjust)
        
       prop_adjust=0; 
    end
alpha=min(1,exp(post_prop-old_posterior+adjustment( posterior_draw,posterior_draw_prop,setup ))*prop_adjust);
elseif bb<setup.number_blocks %a non-dirichlet block
alpha=min(1,exp(post_prop-old_posterior+adjustment( posterior_draw,posterior_draw_prop,setup )));
elseif bb==setup.number_blocks %the weight block
    
    prop_adjust=exp(dirpdf_log(posterior_draw(setup.weight_index),setup.dirichlet_scaling*posterior_draw_prop(min(setup.weight_index):end)/sum(posterior_draw_prop(min(setup.weight_index):end))))/exp(dirpdf_log(posterior_draw_prop(setup.weight_index),setup.dirichlet_scaling*posterior_draw(min(setup.weight_index):end)/sum(posterior_draw(min(setup.weight_index):end))));
   if isnan(prop_adjust)
        
       prop_adjust=0; 
    end
alpha=min(1,exp(post_prop-old_posterior+adjustment( posterior_draw,posterior_draw_prop,setup ))*prop_adjust);
end





if rand<alpha
   posterior_draw=posterior_draw_prop;
   old_posterior=post_prop;
   
   
   if bb==1 || bb==setup.number_blocks
   log_post=log_post_prop;
   
   else
   log_post(bb-1)=post_prop_mm;
   end
   
   acceptances(bb)=acceptances(bb)+1;
   temp_draw=inv_transform( posterior_draw,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );
   xdraw=x_prop;
   add_matrix=add_matrix_prop;
end



   end
if ((nn)/setup.keep_draw)==floor(((nn)/setup.keep_draw))
param_draws(:,(nn)/setup.keep_draw)=temp_draw;
log_posteriors((nn)/setup.keep_draw)=old_posterior;

%later, I'll modify this so we can save the estimated states for the
%different models


%x_draws(:,:,(nn)/setup.keep_draw)=xdraw;



if setup.add_matrices
add_matrices_draws(:,:,(nn)/setup.keep_draw)=add_matrix;
end


end

end


acc_rate=acceptances/setup.number_of_draws;

end
