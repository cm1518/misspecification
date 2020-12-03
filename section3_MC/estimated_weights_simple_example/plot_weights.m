%plots prior vs posterior model weights

%first generate prior draws



for jj=1:length(log_posteriors)
   prior_draw(jj,:)= drchrnd(setup.dirichlet_prior_parameters,1); 
    
    
    
end


figure;
subplot(2,1,1)
histogram(prior_draw(:,1),50)
hold on
histogram(draws_ew(end-1,:),50)
hold off
title('weight on AR(1) model, true model MA(2)')
grid on
subplot(2,1,2)
histogram(prior_draw(:,2),50)
hold on
histogram(draws_ew(end,:),50)
hold off
legend('prior','posterior')
title('weight on MA(1) model, true model MA(2)')
grid on