function [ lnp ] = prior( param,setup )
%computes log prior for normal and gamma prior distributions

if ~isempty(setup.index_normal)
log_normal=(log(normpdf(param(setup.index_normal),setup.normal_prior_means,setup.normal_prior_std)));
else
log_normal=[];
end


if ~isempty(setup.index_gamma)
log_gamma=(log(pdf('Gamma',param(setup.index_gamma),setup.gamma_prior_shape,setup.gamma_prior_scale)));
else
log_gamma=[];
end

%if ~isempty(setup.index_uniform)
%log_uniform=(log(pdf('Uniform',param(setup.index_uniform),setup.uniform_prior_shape,setup.uniform_prior_std)));
%else
%log_uniform=[];
%end

lnp=zeros(size(setup.initial_parameter));
lnp(setup.index_normal)=log_normal;
lnp(setup.index_gamma)=log_gamma;
%lnp(setup.index_uniform)=log_uniform;


end

