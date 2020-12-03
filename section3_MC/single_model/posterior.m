function [ lnp] = posterior( params,setup,data )
%function that returns log posterior


prior_values=prior(params,setup);




if setup.likelihood==1
      % ll_function=SSKF_wrap( params,setup,data );
      [lnpost x add_matrices]=SSKF_wrap_withstates( params,setup,data );
       
       
   elseif setup.likelihood==2
       %ll_function=KF_wrap;
       [lnpost x add_matrices]=KF_wrap_withstates( params,setup,data );
       
   
  
   end

   prior_sum=sum(prior_values);


   lnp=prior_sum+lnpost;
   
end

