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

log_dir=(dirpdf_log(param(setup.weight_index),setup.dirichlet_prior_parameters));
if isnan(log_dir)
    log_dir=-inf;
end


lnp=zeros(size(setup.initial_parameter));
lnp(setup.index_normal)=log_normal;
lnp(setup.index_gamma)=log_gamma;
lnp(end)=log_dir;

end

