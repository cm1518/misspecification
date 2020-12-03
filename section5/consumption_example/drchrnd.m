function r = drchrnd(a,n)
% take a sample from a dirichlet distribution
%a has to be a column vector
p = length(a);
r = gamrnd(repmat(a,1,n),1,p,n);
r = r ./ repmat(sum(r,1),p,1);