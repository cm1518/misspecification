function [ llk ] = KF_wrap( params,setup,data )
[ param ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );   

if setup.initial_provided==1
wrapper_func=str2func(setup.wrapper);
[A B C D]=wrapper_func(param,setup,data);

if setup.TVKF==1
    [ llk xest] = TVKF(A,B,C,D,setup.state_initial,setup.cov_initial,data);
else
[ llk xest] = KF(A,B,C,D,setup.state_initial,setup.cov_initial,data);
end

else
wrapper_func=str2func(setup.wrapper);
[A B C D initial_state initial_cov]=wrapper_func(param,setup,data);

if setup.TVKF==1
[ llk xest] = TVKF(A,B,C,D,initial_state,initial_cov,data);
else
[ llk xest] = KF(A,B,C,D,initial_state,initial_cov,data);
end
    
end
end

