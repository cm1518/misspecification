function [A, B, C_final, D_final, add_matrices] = liquid( params,setup,data )

rho=params(1);
sigma_e1=params(2);
sigma_e2=params(3);
sigma_e3=params(4);
gamma=params(5);    
p_bar=params(6);
d1_bar=params(7);
c1_bar=params(8);
a1_bar=params(9);
d2_bar=params(10);
c2_bar=params(11);
a2_bar=params(12);


%calibrated parameters
omega=0.2;
delta=0.025;
phi=0.95;
beta_1=0.98; %setting to 0.99 makes beta_2>1
r=1/beta_1-1;
beta_2=beta_1*1.02;
xi_1=(1-(1-delta)/(1+r));
xi_2=(1-beta_2*(1-delta-phi*(1+r))-1);



%variables
y=1;
c=2;
d=3;
a=4;
c_1h=5;
c_2h=6;
a_1h=7;
a_2h=8;
d_1h=9;
d_2h=10;
p=11;
yT=12;
yP=13;
c_1ht=14;
c_2ht=15;
a_1ht=16;
a_2ht=17;
d_1ht=18;
d_2ht=19;
Ec_1ht=20;
Ec_2ht=21;
Ed_1ht=22;
Ed_2ht=23;
Ep=24;
mu=25;

g0=zeros(25,25);
g1=zeros(25,25);
C=zeros(25,1);
pi_matrix=zeros(25,5);
PSI=zeros(25,3);


%eq 116
 
g0(1,c_1ht)=1;
g0(1,d_1ht)=-1;
g0(1,p)=-1/xi_1;
g0(1,Ep)=1/xi_1*(1-delta)/(1+r);

%eq 117
g0(2,a_1ht)=1;
g0(2,d_1ht)=-1;
g0(2,p)=-1;

% eq 118
g0(3,a_1ht)=a1_bar;
g0(3,yT)=-1;
g0(3,c_1ht)=c1_bar;
g0(3,p)=p_bar*d1_bar*delta;
g0(3,d_1ht)=p_bar*d1_bar*delta;
g1(3,a_1ht)=a1_bar*(1+r);
PSI(3,2)=-((1+r)*a1_bar-(1-delta)*p_bar*d1_bar)*sigma_e2;


%eq 119
 
g0(4,c_2ht)=1;
g0(4,d_2ht)=-1;
g0(4,p)=-1/xi_2*(1+phi*(beta_2*(1+r)-1));
g0(4,Ep)=beta_2*(1-delta);
g0(4,c_2ht)=(gamma-1)*beta_2*(phi*(1+r)-(1-delta))*c2_bar;
g0(4,d_2ht)=-(gamma-1)*beta_2*(phi*(1+r)-(1-delta))*d2_bar;
g0(4,Ec_2ht)=-(gamma-1)*beta_2*(phi*(1+r)-(1-delta))*c2_bar;
g0(4,Ed_2ht)=(gamma-1)*beta_2*(phi*(1+r)-(1-delta))*d2_bar;


% eq 120
g0(5,p)=p_bar*d2_bar*delta;
g0(5,d_2ht)=p_bar*d2_bar*delta;
g0(5,c_2ht)=c2_bar;
g0(5,a_2ht)=a2_bar;
g0(5,yT)=-1;
g1(5,a_2ht)=(1+r)*a2_bar;

PSI(5,2)=-((1+r)*a2_bar-(1-delta)*p_bar*d2_bar)*sigma_e2;


%   eq 121
g0(6,a_2ht)=a2_bar;
g0(6,p)=1-a2_bar;   
g0(6,d_2ht)=1-a2_bar;

%eq 123

g0(7,c_1h)=1;
g0(7,c_1ht)=-1;
g0(7,yP)=-1;


g0(8,c_2h)=1;
g0(8,c_2ht)=-1;
g0(8,yP)=-1;

%eq 124

g0(9,a_1h)=1;
g0(9,a_1ht)=-1;
g0(9,yP)=-1;




%eq 125

g0(10,d_1h)=1;
g0(10,d_1ht)=-1;
g0(10,yP)=-1;


g0(11,d_2h)=1;
g0(11,d_2ht)=-1;
g0(11,yP)=-1;

%eq 126

g0(12,yP)=1;
g1(12,yP)=1;
PSI(12,2)=sigma_e2;

%eq 127
g0(13,yT)=1;
g1(13,yT)=rho;
PSI(13,1)=sigma_e1;


%eq 128

g0(14,p)=1;
g1(14,p)=phi;
PSI(14,3)=sigma_e3;

%eq 129

g0(15,y)=1;
g0(15,yP)=-1;
g0(15,yT)=-1;

%eq 130

g0(16,c)=1;
g0(16,c_1h)=-omega;
g0(16,c_2h)=-(1-omega);

%eq 131

g0(17,a)=1;
g0(17,a_1h)=-omega;
g0(17,a_2h)=-(1-omega);

%eq 132

g0(18,d)=1;
g0(18,d_1h)=-omega;
g0(18,d_2h)=-(1-omega);

%expectations


g0(19,c_1ht)=1;
g1(19,Ec_1ht)=1;
pi_matrix(19,1)=1;

g0(20,c_2ht)=1;
g1(20,Ec_2ht)=1;
pi_matrix(20,2)=1;


g0(21,d_1ht)=1;
g1(21,Ed_1ht)=1;
pi_matrix(21,3)=1;

g0(22,d_2ht)=1;
g1(22,Ed_2ht)=1;
pi_matrix(22,4)=1;

%eq. 122
g0(23,c_2ht)=(1-1/(beta_2*(1+r)-1)*(beta_2*(1+r)))*(gamma-1)*c2_bar;
g0(23,d_2ht)=(-1+1/(beta_2*(1+r)-1)*(beta_2*(1+r)))*(gamma-1)*d2_bar;
g0(23,mu)=-1;
g0(23,Ec_2ht)=(1/(beta_2*(1+r)-1)*(beta_2*(1+r)))*(gamma-1)*c2_bar;
g0(23,Ed_2ht)=(-1/(beta_2*(1+r)-1)*(beta_2*(1+r)))*(gamma-1)*d2_bar;

%Ep
g0(24,p)=1;
g1(24,Ep)=1;
pi_matrix(24,4)=1;


g0(25,a_2h)=1;
g0(25,a_2ht)=-1;
g0(25,yP)=-1;

%call gensys
 
 [C,Constant,D_sqr,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,C,PSI,pi_matrix);
 D=D_sqr*D_sqr';
if eu(1)~=1 || eu(2)~=1
    C=zeros(size(C)); 
    D=eps+zeros(size(C)); 
 end
  %now add intercept term
 C_final=zeros(26,26);
 C_final(1:25,1:25)=C;
 C_final(26,1:25)=Constant;
 C_final(26,26)=1;
 D_final=zeros(size(D,1)+1,size(D,2)+1);
 D_final(1:size(D,1),1:size(D,2))=D;
 
 
 
 
 A=zeros(3,26); 
 A(1,y)=1;  
 A(2,a)=1;
 A(3,c)=1;


 
 B=zeros(3,1);



B=B*B';

add_matrices=[];