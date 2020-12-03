function [ acc_rate ,last_param_draw,new_posterior] = acceptance_rate( posterior_draw,cov_matrix,scaling,number_of_draws,old_posterior,setup, data,block )
%calculates the acceptance rate of running the standard MH algorithm for
%number_of_draws
acceptances=0;


   
  

for nn=1:number_of_draws
  
param_prop=proposal_draw(posterior_draw(setup.index_block{block}),cov_matrix, scaling,setup);
param_prop=param_prop';
posterior_draw_prop=posterior_draw;
posterior_draw_prop(setup.index_block{block})=param_prop;

post_prop=posterior(posterior_draw_prop,setup,data);
alpha=min(1,exp(post_prop-old_posterior+adjustment( posterior_draw,posterior_draw_prop,setup )));
if rand<alpha
   posterior_draw=posterior_draw_prop;
   old_posterior=post_prop;
   acceptances=acceptances+1;
end

end


acc_rate=acceptances/number_of_draws;
last_param_draw=posterior_draw;
new_posterior=old_posterior;
end

