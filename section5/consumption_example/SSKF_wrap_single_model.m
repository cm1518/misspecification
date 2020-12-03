function [ llk ] = SSKF_wrap_single_model( params,setup,data,mm )

if setup.initial_provided==1


wrapper_func=str2func(setup.wrapper{mm});
[A B C D add_matrices]=wrapper_func([params(setup.index_model{1});params(setup.index_model{mm+1})],setup,data{mm});
[ llk xest] = SSKF(A,B,C,D,setup.state_initial{mm},data{mm});


else

wrapper_func=str2func(setup.wrapper{mm});
[A B C D initial_state initial_cov add_matrices]=wrapper_func([params(setup.index_model{1});params(setup.index_model{mm+1})],setup,data{mm});;
[ llk xest] = SSKF(A,B,C,D,initial_state,data{mm});

    





end

xest=[];
add_matrices=[]; %keep those empty for now. Change this later if we want to store estimates states or something else
