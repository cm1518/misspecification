

n=100000;

draws_prior=drchrnd(setup.dirichlet_prior_parameters,n);
Xi=0:.01:1;
models={'AR(1)','MA(1)'}
figure;
for jj=1:setup.number_models

[F_post]=ksdensity(draws(end-setup.number_models+jj,:),Xi);
[F_prior]=ksdensity(draws_prior(jj,:),Xi);
subplot(setup.number_models,1,jj)

plot(Xi, F_post,Xi,F_prior,'LineWidth',2)
legend('Posterior','Prior')
grid on
title(models{jj})
end

print -depsc