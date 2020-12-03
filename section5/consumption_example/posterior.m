function [ lnp] = posterior( params,setup,data )
%function that returns log posterior
[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );


prior_values=prior(params,setup);


model_weights=params(end-setup.number_models+1:end);



if setup.likelihood==1
      % ll_function=SSKF_wrap( params,setup,data );
      [lnpost x add_matrices]=SSKF_wrap_withstates( params,setup,data );
       
       
   elseif setup.likelihood==2
       %ll_function=KF_wrap;
       [lnpost x add_matrices]=KF_wrap_withstates( params,setup,data );
       
   
  
   end

   prior_sum=sum(prior_values(setup.index_model{1}));
ll_value=0;
for mm=1:setup.number_models
prior_sum=prior_sum+sum(prior_values(setup.index_model{mm+1}))*model_weights(mm);
ll_value=ll_value+lnpost(mm)*model_weights(mm);
end
prior_sum=prior_sum+prior_values(end); %adding the prior for the weights;
lnp=prior_sum+ll_value;
  lnp=lnp-prior_adjustment( params,setup ); %making sure all priors integrate to 1 
end

