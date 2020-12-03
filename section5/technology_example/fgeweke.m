function l=fgeweke(x,mu,prec,chi,A)

% Weigthing function proposed by Geweke (1998) to compute 
% marginal likelihoods.
% This function check the truncation and evaluates the kernel
% of the normal multivariate pdf

dev=(x-mu)*prec*((x-mu)');

if (dev<chi) 
   l=A*exp(-0.5*dev);
else
   l=0;
end
