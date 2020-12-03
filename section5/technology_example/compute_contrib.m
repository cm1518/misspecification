function contrib_array = compute_contrib(ind_vars,total_var)

contrib_array=NaN(size(ind_vars,1),size(ind_vars,3));

for ss=1:size(ind_vars,3)
    contrib_array(:,ss)=diag(ind_vars(:,:,ss))./total_var;
end

end

