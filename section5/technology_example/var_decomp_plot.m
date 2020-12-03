clear;
close all;
clc;
number_draws=500;

horizon=60;






load results_technology_QE
for dd=1:number_draws

ind=randi(size(draws,2),1);
post_med_CL=draws(:,ind);

post_med_JPT_CL=post_med_CL([1:2 14:46]);
[A, B ,C_final, D_final ,add_matrices,inv_shock_std] = modelJPT_for_irf( post_med_JPT_CL,[],[] );
[var_decomp_array_JPT2(:,:,:,dd)] = var_decomp(C_final,D_final,horizon);

end

shocks=[1 2];


variables=[10 58];

%load corresponding number from single model estimation
load JPT_ML;



JPT_lower2=prctile(var_decomp_array_JPT2,5,4);
JPT_upper2=prctile(var_decomp_array_JPT2,95,4);
JPT_med2=prctile(var_decomp_array_JPT2,50,4);




figure;
plot(1:horizon,squeeze(JPT_med(variables(2),shocks(2),:)),'b',1:horizon,squeeze(JPT_med2(variables(2),shocks(2),:)),'r','LineWidth',2)
hold on
plot(1:horizon,squeeze(JPT_lower(variables(2),shocks(2),:)),'b--',1:horizon,squeeze(JPT_lower2(variables(2),shocks(2),:)),'r--','LineWidth',2)
hold on
plot(1:horizon,squeeze(JPT_upper(variables(2),shocks(2),:)),'b--',1:horizon,squeeze(JPT_upper2(variables(2),shocks(2),:)),'r--','LineWidth',2)

grid on
legend('individual model estimation','CL','Location','best')

print -depsc
savefig('var_decomp_QE')

