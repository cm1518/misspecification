function [ acc_rate ,last_param_draw,new_posterior] = acceptance_rate( posterior_draw,cov_matrix,scaling,number_of_draws,old_posterior,setup, data,block )
%calculates the acceptance rate of running the standard MH algorithm for
%number_of_draws
%either all parameters are estimated in the same block or the weights (And only the weights) are in the last block
acceptances=0;


if block ~=1 && block ~=setup.number_blocks
[log_post, ~, ~]=posteriorwithstates_single_model(posterior_draw,setup,data,block-1);

end 
  

for nn=1:number_of_draws
  
  
  if length(setup.index_block{block})==setup.length_param_vector %there is only one block
  param_prop=proposal_draw(posterior_draw(1:(min(setup.weight_index)-1)),cov_matrix, scaling,setup); %random walk proposal
  param_prop=param_prop';
  param_prop_dirichlet=drchrnd(setup.dirichlet_scaling*posterior_draw((min(setup.weight_index):end))/sum(posterior_draw((min(setup.weight_index):end))),1);
  
  param_prop=[param_prop;param_prop_dirichlet];
  elseif block<setup.number_blocks %a non-dirichlet block
  param_prop=proposal_draw(posterior_draw(setup.index_block{block}),cov_matrix, scaling,setup);
  param_prop=param_prop';
  elseif block==setup.number_blocks %the weight block
  param_prop_dirichlet=drchrnd(setup.dirichlet_scaling*posterior_draw((min(setup.weight_index):end))/sum(posterior_draw((min(setup.weight_index):end))),1);
  
  param_prop=param_prop_dirichlet;
  end

posterior_draw_prop=posterior_draw;
posterior_draw_prop(setup.index_block{block})=param_prop;


if block==1 || block==setup.number_blocks;

[post_prop, ~, ~,~]=posteriorwithstates(posterior_draw_prop,setup,data);

else
[post_prop_mm, x_prop, add_matrix_prop]=posteriorwithstates_single_model(posterior_draw_prop,setup,data,block-1);
%post_prop=posterior(param_prop,setup,data);
post_prop=old_posterior-log_post+post_prop_mm;
% if block==2
%    post_prop-posterior(posterior_draw_prop,setup,data) 
% end
end




if length(setup.index_block{block})==setup.length_param_vector
prop_adjust=exp(dirpdf_log(posterior_draw(setup.weight_index),setup.dirichlet_scaling*posterior_draw_prop(min(setup.weight_index):end)/sum(posterior_draw_prop(min(setup.weight_index):end))))/exp(dirpdf_log(posterior_draw_prop(setup.weight_index),setup.dirichlet_scaling*posterior_draw(min(setup.weight_index):end)/sum(posterior_draw(min(setup.weight_index):end))));

alpha=min(1,exp(post_prop-old_posterior+adjustment( posterior_draw,posterior_draw_prop,setup ))*prop_adjust);
elseif block<setup.number_blocks %a non-dirichlet block
alpha=min(1,exp(post_prop-old_posterior+adjustment( posterior_draw,posterior_draw_prop,setup )));
elseif block==setup.number_blocks %the weight block
   prop_adjust=exp(dirpdf_log(posterior_draw(setup.weight_index),setup.dirichlet_scaling*posterior_draw_prop(min(setup.weight_index):end)/sum(posterior_draw_prop(min(setup.weight_index):end))))/exp(dirpdf_log(posterior_draw_prop(setup.weight_index),setup.dirichlet_scaling*posterior_draw(min(setup.weight_index):end)/sum(posterior_draw(min(setup.weight_index):end))));
  if isnan(prop_adjust)
        
       prop_adjust=0; 
   end
alpha=min(1,exp(post_prop-old_posterior+adjustment( posterior_draw,posterior_draw_prop,setup ))*prop_adjust);
% alpha
% post_prop-old_posterior
% posterior_draw_prop(setup.weight_index)
% posterior_draw(setup.weight_index)

end
if rand<alpha
   posterior_draw=posterior_draw_prop;
   old_posterior=post_prop;
   acceptances=acceptances+1;
   if block~=1 && block~=setup.number_blocks
   log_post=post_prop_mm;
   end
end

end


acc_rate=acceptances/number_of_draws;
last_param_draw=posterior_draw;
new_posterior=old_posterior;
end

