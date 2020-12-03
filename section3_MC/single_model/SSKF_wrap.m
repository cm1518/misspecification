function [ llk ] = SSKF_wrap( params,setup,data )
[ param ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );   

if setup.initial_provided==1
wrapper_func=str2func(setup.wrapper);
[A B C D]=wrapper_func(param,setup,data);
[ llk xest] = SSKF(A,B,C,D,setup.state_initial,data);


else
wrapper_func=str2func(setup.wrapper);
[A B C D initial_state initial_cov]=wrapper_func(param,setup,data);
[ llk xest] = SSKF(A,B,C,D,initial_state,data);

    
end
end