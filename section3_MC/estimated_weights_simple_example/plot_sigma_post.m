load results_flat_prior

[F_flat Xi_flat]=ksdensity(draws(1,:));
MSE_flat=mean((draws(1,:)-1).^2);


load results_original_prior
[F_org Xi_org]=ksdensity(draws(1,:));
MSE_org=mean((draws(1,:)-1).^2);

figure;
plot(Xi_flat,F_flat,Xi_org,F_org,'LineWidth',2)
grid on
legend('flat prior','non-flat prior')
title('Posterior of \sigma')
print -depsc

