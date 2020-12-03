function [ lnp x add_matrices] = posteriorwithstates( params,setup,data )
%function that returns log posterior and unobserved states
[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );


prior_values=prior(params,setup);




if setup.likelihood==1
      % ll_function=SSKF_wrap( params,setup,data );
      [lnpost x add_matrices]=SSKF_wrap_withstates( params,setup,data );
       
       
   elseif setup.likelihood==2
       %ll_function=KF_wrap;
       [lnpost x add_matrices]=KF_wrap_withstates( params,setup,data );
       
   
  
   end

   prior_sum_common=sum(prior_values(setup.index_model{1}));
   prior_sum=prior_sum_common;
ll_value=0;
for mm=1:setup.number_models
prior_sum=prior_sum+sum(prior_values(setup.index_model{mm+1}))*setup.model_weights(mm);
ll_value=ll_value+lnpost(mm)*setup.model_weights(mm);
end

   lnp=prior_sum+ll_value;
 
add_matrices=[];
for mm=1:setup.number_models
add_matrices=[add_matrices prior_sum_common+sum(prior_values(setup.index_model{mm+1}))+lnpost(mm)];
end


end


