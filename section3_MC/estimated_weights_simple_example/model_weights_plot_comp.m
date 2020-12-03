load results_equal_prior
draws_equal=draws;
load results_unequal_prior

Xi=0:.01:1;
models={'Small NK, no wages','Small NK','Medium-size NK','Matching'}
figure;
for jj=1:setup.number_models

[F_post1]=ksdensity(draws_equal(end-setup.number_models+jj,:),Xi);
[F_post2]=ksdensity(draws(end-setup.number_models+jj,:),Xi);
subplot(setup.number_models,1,jj)

plot(Xi, F_post1,Xi, F_post2,'LineWidth',2)
legend('Flat prior','Non-Flat Prior')
grid on
title(models{jj})
end

print -depsc