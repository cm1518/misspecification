function [ pa ] = prior_adjustment( param,setup )

%computes (log of) normalizing constant so that priors for restricted variables integrate to 1


%initialization
pa=0; 

if setup.length_log>0 && length(setup.index_normal)>0
[normal_log,index_normal,index_log]=intersect(setup.index_normal,setup.index_log);
normal_log_means=setup.normal_prior_means(index_normal);
normal_log_std=setup.normal_prior_std(index_normal);

normal_log_pa=sum(log(1-normcdf(0,normal_log_means,normal_log_std)));

pa=pa+normal_log_pa;

end


if setup.length_logit>0 && length(setup.index_normal)>0
[normal_logit,index_normal,index_logit]=intersect(setup.index_normal,setup.index_logit);
normal_logit_means=setup.normal_prior_means(index_normal);
normal_logit_std=setup.normal_prior_std(index_normal);

normal_logit_pa=sum(log(normcdf(1,normal_logit_means,normal_logit_std)-normcdf(0,normal_logit_means,normal_logit_std)));

pa=pa+normal_logit_pa;

end


if setup.length_logit_general>0 && length(setup.index_normal)>0
[normal_logit_general,index_normal,index_logit]=intersect(setup.index_normal,setup.index_logit_general);
normal_logit_general_means=setup.normal_prior_means(index_normal);
normal_logit_general_std=setup.normal_prior_std(index_normal);

normal_logit_general_pa=sum(log(normcdf(setup.logit_general_ub(index_logit),normal_logit_general_means,normal_logit_general_std)-normcdf(setup.logit_general_lb(index_logit),normal_logit_general_means,normal_logit_general_std)));

pa=pa+normal_logit_general_pa;

end



if setup.length_logit>0 && length(setup.index_gamma)>0
[gamma_logit,index_gamma,index_logit]=intersect(setup.index_gamma,setup.index_logit);
gamma_logit_shape=setup.gamma_prior_shape(index_gamma);
gamma_logit_scale=setup.gamma_prior_scale(index_gamma);

gamma_logit_pa=sum(log(gamcdf(1,gamma_logit_shape,gamma_logit_scale)));

pa=pa+gamma_logit_pa;

end

if setup.length_logit_general>0 && length(setup.index_gamma)>0
[gamma_logit_general,index_gamma,index_logit]=intersect(setup.index_gamma,setup.index_logit_general);
gamma_logit_general_shape=setup.gamma_prior_shape(index_gamma);
gamma_logit_general_scale=setup.gamma_prior_scale(index_gamma);

gamma_logit_general_pa=sum(log(gamcdf(setup.logit_general_ub(index_logit),gamma_logit_general_shape,gamma_logit_general_scale)-gamcdf(setup.logit_general_ub(index_logit),gamma_logit_general_shape,gamma_logit_general_scale)));

pa=pa+gamma_logit_general_pa;

end
