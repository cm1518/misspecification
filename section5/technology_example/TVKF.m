function [ llk xest] = TVKF(A,B,C,D,x10,S10,Y)
%time varying Kalman Filter, returns estimates of the state and the log-likelihood
%follows Hamilton, pages 380-381 (in particular, equations
%(13.2.19),(13.2.20),(13.2.22) and(13.4.1)).
%state space system is given by:
%y_t=A_tx_t+u^y_t
%x_t+1=C_tx_t+u^x_t
%u^y_t N(0,B_t)
%u^x_t N(0,D_t)
%We assume that the observables Y are of dimension N by T
%initial state is x10 (x(t|t-1)). Same notation for the covariance matrix
%of the state
[N T]=size(Y);
[~, sizex]=size(A);
xest=zeros(sizex,T+1);
xest(:,1)=x10;
llk=0;
for i=1:T
    
    A_temp=A(:,:,i);
    B_temp=B(:,:,i);
    C_temp=C(:,:,i);
    D_temp=D(:,:,i);
    
    temp=A_temp*S10*A_temp'+B_temp;
    temp2=Y(:,i)-A_temp*xest(:,i);
    llk=llk+log(2*pi)*(-N/2)+logdet(temp, 'chol')*(-.5)-.5*(temp2)'/temp*temp2;
    K=C_temp*S10*A_temp'/(temp);
    xest(:,i+1)=C_temp*xest(:,i)+K*temp2;
    S10=C_temp*(S10-S10*A_temp'/temp*A_temp*S10)*C_temp'+D_temp;
    S10=(S10+S10')/2; %making sure the covariance matrix stays symmetric
end

