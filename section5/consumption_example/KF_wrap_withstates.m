function [ llk xest add_matrices] = KF_wrap_withstates( params,setup,data )
%add_matrices ordered last in wrapper file (return empty matrix if not
%needed)

llk=zeros(setup.number_models,1);
if setup.initial_provided==1

for mm=1:setup.number_models
 try
wrapper_func=str2func(setup.wrapper{mm});
[A B C D add_matrices]=wrapper_func([params(setup.index_model{1});params(setup.index_model{mm+1})],setup,data{mm});

if setup.TVKF==1
    [ llk(mm) xest] = TVKF(A,B,C,D,setup.state_initial{mm},setup.cov_initial{mm},data{mm});
else
[ llk(mm) xest] = KF(A,B,C,D,setup.state_initial{mm},setup.cov_initial{mm},data{mm});
end
 catch
     llk(mm)=-inf;
 end

end

else
for mm=1:setup.number_models
    try
wrapper_func=str2func(setup.wrapper{mm});
[A B C D initial_state initial_cov add_matrices]=wrapper_func([params(setup.index_model{1});params(setup.index_model{mm+1})],setup,data{mm});;
if setup.TVKF==1
[ llk(mm) xest] = TVKF(A,B,C,D,initial_state ,initial_cov,data{mm});
else
[ llk(mm) xest] = KF(A,B,C,D,initial_state,initial_cov,data{mm});
end
    catch
        llk(mm)=-inf;
    end
end



end

xest=[];
add_matrices=[]; %keep those empty for now. Change this later if we want to store estimates states or something else
