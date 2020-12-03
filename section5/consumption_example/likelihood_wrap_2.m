function [lnpost]=likelihood_wrap_2( params,setup,data )

[ params ] = [setup.alt_start_value(1:16);params(1:setup.size_obs);setup.alt_start_value(17:end);params(setup.size_obs+1:end)];


try
[ lnpost, ~, ~, ~] = likelihood(data, params,setup );
catch ME
    
   lnpost=-1e100;
end
    
if isnan(lnpost)
    lnpost=-1e100;
end


end

