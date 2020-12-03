function [A B C_final D_final add_matrices] = nk1( params,setup,data )
%code solves the NK model using gensys and observables x (output gap), pi, and r
%model now features three structural shocks

%unwrap parameters
rho_u=params(1);
phi=params(2);
psix=params(3);
psir=params(4);
psip=params(5);
sigma_e=params(6);
sigma_u=params(7);
sigma_g=params(8);
rho_g=params(9);



pi_bar=params(11);
theta=params(10);

y_bar=1.67;

%fixing other coefficients (values from Fabio's code)
betta= 0.99;
r_bar=1/betta;

%setting up system matrices
% variables are ordered as [x p r E_t(x_(t+1)) E_t(p_(t+1)) u g x(-1)]

g0=zeros(8,8);
g0(1,1)=1;
g0(1,3)=1/phi;
g0(1,4)=-1;
g0(1,5)=-(1/phi);
g0(1,6)=-(1-rho_u);
g0(1,7)=-(1/phi)*rho_g;
g0(2,2)=1;
g0(2,5)=-betta/(1+betta);
g0(2,1)=-theta/(1+betta);
g0(2,6)=theta/(1+betta);
g0(3,3)=1;
g0(3,1)=-(1-psir)*psix;
g0(3,2)=-(1-psir)*psip;
g0(4,1)=1;
g0(5,2)=1;
g0(6,6)=1;
g0(7,7)=1;
g0(8,8)=1;

g1=zeros(8,8);
g1(2,2)=1/(1+betta);
g1(3,3)=psir;
g1(4,4)=1;
g1(5,5)=1;
g1(6,6)=rho_u;
g1(7,7)=rho_g;
g1(8,1)=1;

c=zeros(8,1);

pi_matrix=zeros(8,2);
pi_matrix(4,1)=1;
pi_matrix(5,2)=1;

psi=zeros(8,3);
psi(3,1)=sigma_e;
psi(6,2)=sigma_u;
psi(7,3)=sigma_g;

 %call gensys
 
 [C,Constant,D_sqr,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi_matrix);
 D=D_sqr*D_sqr';

 %now add intercept term
 C_final=zeros(9,9);
 C_final(1:8,1:8)=C;
 C_final(9,9)=1;
 D_final=zeros(size(D,1)+1,size(D,2)+1);
 D_final(1:size(D,1),1:size(D,2))=D;
 
 
 
 
 A=zeros(3,9);
 A(1,2)=400; %inflation
 A(1,9)=pi_bar;
 A(2,3)=400;
A(2,9)=(pi_bar+400*(r_bar-1)); %nominal interest rate
 A(3,1)=400;% output growth
 A(3,7)=400;
 A(3,8)=-400;% output growth
 A(3,9)=y_bar; % NB: SS values already in annualized percent
 B=zeros(3,3);



add_matrices=[];