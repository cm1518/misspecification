function [A B C_final D_final add_matrices] = RRR( params,setup,data )
%code solves the model with sticky prices and wages as well as price indexation from the JME (2005) paper
%by Rabanal and Rubio-Ramirez. All data are detrended before estimation, hence there are no constants floating around.

%unwrap parameters
sigma=params(1);
gamma=params(2);
rho_r=params(3);
gamma_pi=params(4);
gamma_y=params(5);
rho_a=params(6);
sigma_a=params(7);
rho_g=params(8);
sigma_g=params(9);
sigma_z=params(10);
sigma_lambda=params(11);
gamma_b=params(12);

kappa_w=params(13);
kappa_p=params(14);
%fixing other coefficients 
betta= 0.99;
delta=0.36;
gamma_f=betta*(1-betta*gamma_b);
%setting up system matrices
% variables are ordered as [y r pi n mrs w p delta(w) mc g a z lambda E_tpi_t+1 E_t delta(w)_t+1 E_tw_t+1 w(-1) p(-1)]

%indices for positions of variables
y=1;
r=2;
pi=3;
n=4;
mrs=5;
w=6;
p=7;
delta_w=8;
mc=9;
g=10;
a=11;
z=12;
lambda=13;
Epi=14;
Edeltaw=15;
Ey=16;



g0=zeros(16,16);
g1=zeros(16,16);
c=zeros(16,1);
pi_matrix=zeros(16,3);
psi=zeros(16,4);



%Euler equation

g0(1,y)=1;
g0(1,Ey)=-1;
g0(1,r)=sigma;
g0(1,Epi)=-sigma;
g0(1,g)=-(1-rho_g);


%production

g0(2,y)=1;
g0(2,a)=-1;
g0(2,n)=-(1-delta);

%marginal cost


g0(3,mc)=1;
g0(3,w)=-1;
g0(3,p)=1;
g0(3,n)=-1;
g0(3,y)=1;


%mrs

g0(4,mrs)=1;
g0(4,y)=-(1/sigma);
g0(4,n)=-gamma;
g0(4,g)=1;


%taylor rule

g0(5,r)=1;
g1(5,r)=rho_r;
g0(5,y)=-(1-rho_r)*gamma_y;
g0(5,pi)=-(1-rho_r)*gamma_pi;
g0(5,z)=-1;


%real wages

g0(6,w)=1;
g0(6,p)=-1;
g0(6,pi)=1;
g0(6,delta_w)=-1;
g1(6,w)=1;
g1(6,p)=-1;

%inflation
g0(7,pi)=1;
g0(7,p)=-1;
g1(7,p)=-1;



%technology shock

g0(8,a)=1;
g1(8,a)=rho_a;
psi(8,1)=sigma_a;

%g shock

g0(9,g)=1;
g1(9,g)=rho_g;
psi(9,2)=sigma_g;

%z shock

g0(10,z)=1;
psi(10,3)=sigma_z;


%lambda shock

g0(11,lambda)=1;
psi(11,4)=sigma_lambda;

%wage phillips curve

g0(12,delta_w)=1;
g0(12,Edeltaw)=-betta;
g0(12,mrs)=-kappa_w;
g0(12,w)=kappa_w;
g0(12,p)=-kappa_w;

%price phillips curve
g0(13,pi)=1;
g0(13,Epi)=-gamma_f;
g0(13,mc)=-kappa_p;
g0(13,lambda)=-kappa_p;
g1(13,pi)=gamma_b;

%inflation expectations

g0(14,pi)=1;
g1(14,Epi)=1;
pi_matrix(14,1)=1;


%wage growth expectations

g0(15,delta_w)=1;
g1(15,Edeltaw)=1;
pi_matrix(15,2)=1;


%output gap expectations

g0(16,y)=1;
g1(16,Ey)=1;
pi_matrix(16,3)=1;


 %call gensys
 
 [C,Constant,D_sqr,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi_matrix);
 D=D_sqr*D_sqr';

 %now add intercept term
 C_final=zeros(17,17);
 C_final(1:16,1:16)=C;
 C_final(17,17)=1;
 D_final=zeros(size(D,1)+1,size(D,2)+1);
 D_final(1:size(D,1),1:size(D,2))=D;
 
 
 
 
 A=zeros(4,17);
 A(1,pi)=100; %inflation
 
 A(2,r)=100;
 A(3,y)=100;% output gap
 A(4,w)=100;%wage
 B=zeros(4,4);



add_matrices=[];