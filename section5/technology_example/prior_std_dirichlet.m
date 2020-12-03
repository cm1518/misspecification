a_0=sum(temp_prior_weights);

var_prior=(temp_prior_weights.*(a_0-temp_prior_weights))./(a_0^2*(a_0+1));
std_prior=sqrt(var_prior);