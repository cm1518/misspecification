function marginal_dens = m_harmonic(x2,logpo2)

% function marginal = marginal_density()
% Computes the marginal density
%
%
% OUTPUTS
%    marginal_dens: log marginal density
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2005-2007)
% Gnu Public License.
%
% Modified by Mu-Chun Wang (04-2008)


%global options_

warning off
options_.npar=size(x2,2);
options_.nruns=size(x2,1);
lpost_mode = -Inf;
MU=mean(x2(1:end,:));
%SIGMA=cov(x2(1:end,:),1);
SIGMA=cov(x2);

lpost_mode = max(lpost_mode,max(logpo2(1:end)));

%disp(' ');
%disp('MH: I''m computing the posterior log marginal density (modified harmonic mean)... ');
detSIGMA = det(SIGMA);
invSIGMA = inv(SIGMA);
marginal = zeros(9,2);
linee = 0;
check_coverage = 1;
increase = 1;
while check_coverage
  for p = 0.1:0.1:0.9;
    %critval = qchisq(p,options_.npar);
    critval = chi2inv(p,options_.npar);
    
    %critval = chi2inv(p,1);
    %ifil = FirstLine;
    ifil = 1;
    tmp = 0;
    %for n = FirstMhFile:TotalNumberOfMhFiles
      %for b=1:nblck
	%load([ MhDirectoryName '/' M_.fname '_mh' int2str(n) '_blck' int2str(b)],'x2','logpo2');
	EndOfFile = size(x2,1);
	for i = ifil:EndOfFile
	  deviation  = (x2(i,:)-MU)*invSIGMA*(x2(i,:)-MU)';
	  if deviation <= critval
	    lftheta = -log(p)-(options_.npar*log(2*pi)+log(detSIGMA)+deviation)/2;
	    tmp = tmp + exp(lftheta - logpo2(i) + lpost_mode);
	  end
	end
      %end
      %ifil = 1;
    %end
    linee = linee + 1;
    warning off all
    marginal(linee,:) = [p, lpost_mode-log(tmp/options_.nruns)];
    warning on all
  end
  %abs((marginal(9,2)-marginal(1,2))/marginal(9,2))
  %increase
  if abs((marginal(9,2)-marginal(1,2))/marginal(9,2)) > 0.01 | any(isinf(marginal(:,2)))
    if increase == 1
      disp('MH: The support of the weighting density function is not large enough...')
      disp('MH: I increase the variance of this distribution.')
      increase = 1.2*increase;
      invSIGMA = inv(SIGMA*increase);
      detSIGMA = det(SIGMA*increase);
      linee    = 0;   
    else
      %disp('MH: Let me try again.')
      increase = 1.2*increase;
      invSIGMA = inv(SIGMA*increase);
      detSIGMA = det(SIGMA*increase);
      linee    = 0;
      if increase > 20
	check_coverage = 0;
	clear invSIGMA detSIGMA increase;
	%disp('MH: There''s probably a problem with the modified harmonic mean estimator.')    
      end
    end  
  else
    check_coverage = 0;
    clear invSIGMA detSIGMA increase;
    disp('MH: Modified harmonic mean estimator, done!')
  end
end

marginal_dens = mean(marginal(:,2));
marginal