clear all
clc


sixmin = @(x) (x-5)^2;


options = optimoptions(@fminunc,'FiniteDifferenceType','central','MaxIterations',5000,'MaxFunctionEvaluations',50000);

problem = createOptimProblem('fminunc', ...
    'x0',10,'objective',sixmin, ...
    'lb',-10,'ub',10,'options',options);

ms = MultiStart;

[xming,fming,flagg,outptg,manyminsg] = run(ms,problem,k);