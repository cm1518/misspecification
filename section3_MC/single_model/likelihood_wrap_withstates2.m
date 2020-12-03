function [lnpost x add_matrices]=likelihood_wrap_withstates2( params,setup,data )
% [ param ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );   
% 
% current_matrices=zeros(1,1,2);
% current_matrices(1,1,1)=param(2);
% current_matrices(1,1,2)=param(3);
% lagged_matrices=zeros(1,1,2,1);
% lagged_matrices(:,:,1,1)=param(4);
% %lagged_matrices(2,1,1,1)=500    ;
% lagged_matrices(:,:,2,1)=param(5);



[ constant current_matrices lagged_matrices ] = wrap_BM( params, setup );

[ lnpost indicator x uvec] = likelihood(data, constant, current_matrices,lagged_matrices,setup.initial_eps, setup.lags );
add_matrices=x;
%x=0;

if isnan(lnpost)
    lnpost=-inf;
end
end

