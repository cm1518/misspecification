function [var_decomp_array] = var_decomp(C,D,max_horizon)
%computes variance decomposition for first-order VAR up to
%max_horizon
var_decomp_array=NaN(size(D,1),size(D,2),max_horizon);
shock_impact_temp=D;
total_variance=zeros(size(D,1),size(D,1));
ind_shock_var=zeros(size(D,1),size(D,1),size(D,2));

for hh=1:max_horizon-1
    total_variance=shock_impact_temp*shock_impact_temp'+total_variance;
    for ee=1:size(D,2)
   select=zeros(size(D,2),size(D,2));
   select(ee,ee)=1;
   ind_shock_var(:,:,ee)=(shock_impact_temp*select*shock_impact_temp')+ind_shock_var(:,:,ee); %contribution of each shock to total variance
  
   end
    var_decomp_array (:,:,hh)=compute_contrib(ind_shock_var,diag(total_variance));
  
  shock_impact_temp=C*shock_impact_temp; 
end

end

