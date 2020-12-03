function [ llk, xest, add_matrices, l_vec] = KF_wrap_withstates_l_vec( params,setup,data )
%add_matrices ordered last in wrapper file (return empty matrix if not
%needed)

%this code only works if all models use the same sample size (not
%necessarily the same observables, though)
T=size(data{1},2)-1;
l_vec=zeros(T,setup.number_models);
llk=zeros(setup.number_models,1);
if setup.initial_provided==1

for mm=1:setup.number_models
 
    
    
    if mm~=3
    try
wrapper_func=str2func(setup.wrapper{mm});
[A, B, C, D, add_matrices]=wrapper_func([params(setup.index_model{1});params(setup.index_model{mm+1})],setup,data{mm});

if setup.TVKF==1
    [ llk(mm), xest, temp_l] = TVKF_l_vec(A,B,C,D,setup.state_initial{mm},setup.cov_initial{mm},data{mm});
    l_vec(1:end,mm)=temp_l(2:end);
else
[ llk(mm), xest ,temp_l] = KF_l_vec(A,B,C,D,setup.state_initial{mm},setup.cov_initial{mm},data{mm});
l_vec(1:end,mm)=temp_l(2:end);
end
 catch
     llk(mm)=-inf;
    end

 
    else
        
       try
wrapper_func=str2func(setup.wrapper{mm});
[A, B, C, D, add_matrices]=wrapper_func([params(setup.index_model{1});params(setup.index_model{mm+1})],setup,data{mm});

if setup.TVKF==1
    [ llk(mm), xest, l_vec(1:end,mm)] = TVKF_l_vec(A,B,C,D,setup.state_initial{mm},setup.cov_initial{mm},data{mm});
else
[ llk(mm), xest ,l_vec(1:end,mm)] = KF_l_vec(A,B,C,D,setup.state_initial{mm},setup.cov_initial{mm},data{mm});
end
 catch
     llk(mm)=-inf;
 end 
        
        
    end
    
    
end

else
for mm=1:setup.number_models
    try
wrapper_func=str2func(setup.wrapper{mm});
[A, B, C, D, initial_state, initial_cov, add_matrices]=wrapper_func([params(setup.index_model{1});params(setup.index_model{mm+1})],setup,data{mm});;
if setup.TVKF==1
[ llk(mm), xest,l_vec(:,mm)] = TVKF_l_vec(A,B,C,D,initial_state ,initial_cov,data{mm});
else
[ llk(mm), xest, l_vec(:,mm)] = KF_l_vec(A,B,C,D,initial_state,initial_cov,data{mm});
end
    catch
        llk(mm)=-inf;
    end
end



end

xest=[];
add_matrices=[]; %keep those empty for now. Change this later if we want to store estimates states or something else
