function [ lnp ] = postblock( params,setup,data,block, params_estimated )
%computes the posterior as a function of the parameters in one block,
%keeping the parameters in all other blocks fixed
params_total=params_estimated;

if setup.number_blocks>1
params_total(setup.index_block{block})=params;
elseif setup.number_blocks==1
  params_total(1:min(setup.weight_index-1))=params;  
end

[ lnp] = posterior( params_total,setup,data );
end

