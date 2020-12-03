function [A, B, C_final, D_final, add_matrices] = exp_utility( params,setup,data )

rho=params(1);
sigma_e1=params(2);
a_bar=params(3); 
y_bar=params(4); 
c_bar=params(5);  
sigma_e2=params(6);
sigma_e3=params(7);
rho2=params(8);
theta_transform=params(9); %1/(theta*c_bar)
sigma_bar=params(10);

G=0;

yT=1;
c=2;
a=3;
a_total=4;
c_total=5;
Ec=6;
yP=7;
sigma=8;


g0=zeros(8,8);
g1=zeros(8,8);
C=zeros(8,1);
pi_matrix=zeros(8,1);
PSI=zeros(8,3);


%eq 29
 
g0(1,c)=1;
g0(1,Ec)=-1;
g0(1,sigma)=theta_transform;
PSI(1,2)=theta_transform;

% eq 30
g0(2,a)=1;
g0(2,yT)=-y_bar;
g0(2,c)=c_bar;
g0(2,sigma)=a_bar*sigma_bar;
g1(2,a)=a_bar;
PSI(2,2)=-a_bar*sigma_e2;



%   eq 31 
g0(3,yP)=1;
g1(3,yP)=1;
PSI(3,2)=sigma_e2;
C(3,1)=G;
   

%eq 32
g0(4,yT)=1;
g1(4,yT)=rho;
PSI(4,1)=sigma_e1;
%Ec

g0(5,c)=1;
g1(5,Ec)=1;
pi_matrix(5,1)=1;

%eq 34
g0(6,c_total)=1;
g0(6,c)=-1;
g0(6,yP)=-1;

%eq 35
g0(7,a_total)=1;
g0(7,a)=-1;
g0(7,yP)=-1;
 
 %eq 33
g0(8,sigma)=1;
g1(8,sigma)=1;
PSI(8,3)=sigma_e3;


%call gensys
 
 [C,Constant,D_sqr,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,C,PSI,pi_matrix);
 
 D=D_sqr*D_sqr';
if eu(1)~=1 || eu(2)~=1
    C=zeros(size(C)); 
    D=eps+zeros(size(C)); 
 end
  %now add intercept term
 C_final=zeros(9,9);
 C_final(1:8,1:8)=C;
 C_final(9,1:8)=Constant;
 C_final(9,9)=1;
 D_final=zeros(size(D,1)+1,size(D,2)+1);
 D_final(1:size(D,1),1:size(D,2))=D;
 
 
 
 
 A=zeros(3,9); 
 A(1,yP)=1; 
 A(1,yT)=1; 
 A(2,a_total)=1;
 A(3,c_total)=1;

 B=zeros(3,1);


B=B*B';

add_matrices=[];