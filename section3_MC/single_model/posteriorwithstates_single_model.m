function [ lnp x add_matrices] = posteriorwithstates_single_model( params,setup,data ,mm)
%function that returns log posterior and unobserved states
[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );


prior_values=prior(params,setup);




if setup.likelihood==1
      % ll_function=SSKF_wrap( params,setup,data );
      [lnpost x add_matrices]=SSKF_wrap_withstates_single_model( params,setup,data,mm );
       
       
   elseif setup.likelihood==2
       %ll_function=KF_wrap;
       [lnpost x add_matrices]=KF_wrap_withstates_single_model( params,setup,data,mm );
       
   
  
   end

  
ll_value=0;

prior_sum=sum(prior_values(setup.index_model{mm+1}))*setup.model_weights(mm);
ll_value=lnpost*setup.model_weights(mm);


   lnp=prior_sum+ll_value;
 
add_matrices=[];


end


