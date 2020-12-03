load ARMA_050605_T50
prd1=prior_draw(:,1);
pod1=draws_ew(4,:);

load ARMA_100605_T50
prd2=prior_draw(:,1);
pod2=draws_ew(4,:);

load ARMA_100308_T50
prd3=prior_draw(:,1);
pod3=draws_ew(4,:);

load ARMA_100902_T50
prd4=prior_draw(:,1);
pod4=draws_ew(4,:);


 figure;
subplot(2,2,1)
histogram(prd1,100); % 'Normalization', 'pdf')
hold on
histogram(pod1,100); %'Normalization', 'pdf')
hold off
title('DGP 1')
grid on
subplot(2,2,2)
histogram(prd2,100); % 'Normalization', 'pdf')
hold on
histogram(pod2,100); %'Normalization', 'pdf')
hold off
title('DGP 2')
grid on
subplot(2,2,3)
histogram(prd3,100); % 'Normalization', 'pdf')
hold on
histogram(pod3,100); %'Normalization', 'pdf')
hold off
title('DGP 3')
grid on
subplot(2,2,4)
histogram(prd4,100); % 'Normalization', 'pdf')
hold on
histogram(pod4,100); %'Normalization', 'pdf')
hold off
title('DGP 4')
grid on

