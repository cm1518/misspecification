function [A, B, C_final, D_final, add_matrices] = rbc( params,setup,data )

rho=params(1);
sigma_e1=params(2);
c_bar=params(3); 
k_bar=params(4); 
k1_bar=params(5);  
sigma_meas_a=params(6);
sigma_e2=params(7);
gamma=params(8);
N1N=params(9);
eta=params(10);
G=0;
%calibrated parameters
alpha=0.4;
r_bar=1/0.99-1;


yT=1;
c=2; %note that c here is c_1 in the notes
a=3;
a_total=4;
c_total=5;
Ec=6;
yP=7;
y_total=8;
N=9;
Y=10;
r=11;
EN=12;
Er=13;

g0=zeros(13,13);
g1=zeros(13,13);
C=zeros(13,1);
pi_matrix=zeros(13,3);
PSI=zeros(13,2);


%eq 46
 
g0(1,c)=gamma;
g0(1,N)=eta+1;
g0(1,Y)=-1;



% eq 47
g0(2,c)=-gamma;
g0(2,Ec)=gamma;
g0(2,Er)=-r_bar/(1+r_bar);
PSI(2,2)=sigma_e2;
%eq 48

g0(3,r)=1;
g0(3,Y)=-alpha/(1+r_bar);
g1(3,a)=-alpha/(1+r_bar);


%eq 69
g0(4,Y)=1;
g0(4,N)=-(1-alpha);
g0(4,yT)=-(1-alpha);
g1(4,a)=alpha;


%eq 70
g0(5,Y)=1;
g0(5,c)=-c_bar;
g0(5,a)=-k_bar;
g1(5,a)=k1_bar;

PSI(5,2)=-k1_bar*sigma_e2;

%   eq 71 
g0(6,yP)=1;
g1(6,yP)=1;
PSI(6,2)=sigma_e2;

   

%eq 72
g0(7,yT)=1;
g1(7,yT)=rho;
PSI(7,1)=sigma_e1;
%Ec

g0(8,c)=1;
g1(8,Ec)=1;
pi_matrix(8,1)=1;

%EN

g0(9,N)=1;
g1(9,EN)=1;
pi_matrix(9,2)=1;

%Er

g0(10,r)=1;
g1(10,Er)=1;
pi_matrix(10,3)=1;


%eq 73
g0(11,c_total)=1;
g0(11,c)=-1;
g0(11,yP)=-1;

%eq 74
g0(12,a_total)=1;
g0(12,a)=-1;
g0(12,yP)=-1;

%eq 75
g0(13,y_total)=1;
g0(13,Y)=-1;
g0(13,yP)=-1;

%call gensys
 
 [C,Constant,D_sqr,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,C,PSI,pi_matrix);
 D=D_sqr*D_sqr';
% if eu(1)~=1 || eu(2)~=1
%     C=zeros(size(C)); 
%     D=eps+zeros(size(C)); 
%  end
  %now add intercept term
 C_final=zeros(14,14);
 C_final(1:13,1:13)=C;
 %C_final(14,1:13)=Constant;
 C_final(14,14)=1;
 D_final=zeros(size(D,1)+1,size(D,2)+1);
 D_final(1:size(D,1),1:size(D,2))=D;
 
 
 
 
 A=zeros(3,14); 
 A(1,y_total)=1;  
 A(2,a_total)=1;
 A(3,c_total)=1;

 B=zeros(3,1);
B(2,1)=sigma_meas_a;


B=B*B';

add_matrices=[];