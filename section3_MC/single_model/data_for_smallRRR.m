test = xlsread('variables_to_detrend.xlsx');
data=test';

%build X matrix of RHS variables
X=[];
for jj=1:3
   % X=[X [0:setup.sample_size-1]'.^(jj-1)];
    X=[X (2*[0:size(data,2)-1]'/(size(data,2)-1)-1).^(jj-1)]; %rescale polynomials to make algorithm more stable
    % see https://www.rssd.esa.int/SP/LISAPATHFINDER/docs/Data_Analysis/DA_Nine/Annex7b.pdf
end



%run regressions to get trend
intercepts=zeros(2,1);
poly_coefficients=zeros(2,2);
for kk=1:2
   params_temp=X\data(kk,:)'; 
    
  
       intercepts(kk)=params_temp(1);
   
     poly_coefficients(kk,:)=params_temp(2:end);  
   
end
trends=[intercepts poly_coefficients]*X';

figure;
for jj=1:2
    subplot(2,1,jj)
plot(1:size(data,2),data(jj,:)',1:size(data,2),trends(jj,:));
title('data and trends')
end
data_detrended=data-trends;

test = xlsread('variables_to_demean.xlsx');
data=test';

%build X matrix of RHS variables
X=[];
for jj=1:1
   % X=[X [0:setup.sample_size-1]'.^(jj-1)];
    X=[X (2*[0:size(data,2)-1]'/(size(data,2)-1)-1).^(jj-1)]; %rescale polynomials to make algorithm more stable
    % see https://www.rssd.esa.int/SP/LISAPATHFINDER/docs/Data_Analysis/DA_Nine/Annex7b.pdf
end



%run regressions to get trend
intercepts=zeros(2,1);
poly_coefficients=zeros(2,2);
for kk=1:2
   params_temp=X\data(kk,:)'; 
    
  
       intercepts(kk)=params_temp(1);
   
     
   
end
trends=[intercepts]*X';
data_demeaned=data-trends;
figure;
for jj=1:2
    subplot(2,1,jj)
plot(1:size(data,2),data(jj,:)',1:size(data,2),trends(jj,:));
title('data and trends')
end



data=[data_demeaned;data_detrended(1,2:end)];

save data data