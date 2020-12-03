    for i = 1:t;
        BBtempor = Bt_postmean(M+1:end,i);
        BBtempor = reshape(BBtempor,M*p,M)';
        ctemp1 = [BBtempor; eye(M*(p-1)) zeros(M*(p-1),M)];
        max_eig(i)=max(abs(eig(ctemp1)));
           
    end
    
    plot(yearlab,max_eig)