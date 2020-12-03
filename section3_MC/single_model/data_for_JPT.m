% we use the data that accompanies the SW AER paper to estimate the JPT model. Note that JPT use the same model variables in estimation, but the 
%data is constructed differently, which, according to their paper, makes a difference

test = xlsread('sw_data.xlsx');
data=test';


save data data