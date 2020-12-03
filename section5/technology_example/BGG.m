function [A B C_final D_final add_matrices] = BGG( params,setup,data )

kappa_p=params(1);
gamma_b=params(2);
eta=params(3); 
alpha=params(4); 
delta=params(5);  
rhoA=params(6); 
rhoG=params(7); 
psi=params(8); 
%K_N=params(9); 
gamma=params(9); 

rho=params(10); 
S=params(11);
sigma_A=params(12);
sigma_G=params(13);
sigma_MP=params(14);



%calibrated parameters
betta=0.99; 
R         = 1/betta      ;
C_Y	    =	0.61        ; 
Ce_Y	  =	0.01        ; 
I_Y 	  =	0.18    	;
G_Y 	  =	0.20        ;
K_N 	  =	2.00        ;
Y_N 	  =	0.28        ;       
X   	  =	1.10        ;
eps       = (1-delta)/((1-delta) + ((alpha/X)*(Y_N/K_N)));
omega     = 0.99        ;
nu        = 0.5         ;
gamma_f=betta*(1-betta*gamma_b);
%indices for positions of variables

y=1; 
c=2; 
i=3;
g=4; 
ce=5;
n=6;
rk=7;
r=8;
q=9;
k=10; 
x=11; 
a=12; 
h=13; 
pi=14; 
rn=15;
Epi=16;
Ec=17; 
Erk=18;
mp_shock=19;

g0=zeros(19,19);
g1=zeros(19,19);
C=zeros(19,1);
pi_matrix=zeros(19,2);
PSI=zeros(19,3);


% y = C_Y*c + I_Y*i + G_Y*g + Ce_Y*ce;            //4.14
 
g0(1,y)=1;
g0(1,c)=-C_Y;
g0(1,i)=-I_Y;
g0(1,g)=-G_Y;
g0(1,ce)=-Ce_Y;


% c = -r + c(+1);                                 //4.15  
g0(2,c)=1;
g0(2,r)=1;
g0(2,Ec)=-1;

 
%    rn = r + pi(+1);                                //nominal int. rate
g0(3,rn)=1;
g0(3,r)=-1;
g0(3,Epi)=-1;


%    ce = n;                                         //4.16
g0(4,ce)=1;
g0(4,n)=-1;

%    rk(+1) - r = -nu*(n -(q + k));                  //4.17
g0(5,Erk)=1;
g0(5,r)=-1;
g0(5,n)=nu;
g0(5,q)=-nu;
g0(5,k)=-nu;

   % rk = (1-eps)*(y - k(-1) - x) + eps*q - q(-1);   //4.18
g0(6,rk)=1;
g0(6,y)=-(1-eps);
g0(6,x)=(1-eps);
g0(6,q)=-eps;
g1(6,k)=-(1-eps);
g1(6,q)=-1;
   
  % q = psi*(i - k(-1));                            //4.19    
g0(7,q)=1;
g0(7,i)=-psi;
g1(7,k)=-psi;  

  
 %y = a + alpha*k(-1) + (1-alpha)*omega*h;        //4.20
 g0(8,y)=1;
 g0(8,a)=-1;
 g0(8,h)=-(1-alpha)*omega;
 g1(8,k)=alpha;
 

 %y - h - x - c = (eta^(-1))*h;                   //4.21
g0(9,y)=1;
g0(9,h)=-1;
g0(9,x)=-1;
g0(9,c)=-1;
g0(9,h)=-1/eta;  


 %PC                  //4.22
%price phillips curve (NB: x inversely related to mc)
g0(10,pi)=1;
g0(10,Epi)=-gamma_f;
g0(10,x)=kappa_p;
g1(10,pi)=gamma_b;


 
 %k = delta*i + (1-delta)*k(-1);                  //4.23
g0(11,k)=1; 
g0(11,i)=-delta;
g1(11,k)=1-delta;

 %n = gamma*R*K_N*(rk - r(-1)) + r(-1) + n(-1);   //4.24
g0(12,n)=1;
g1(12,rk) = gamma*R*K_N;
g1(12,r)=-gamma*R*K_N+1;
g1(12,n)=1;

 % rn = rho*rn(-1) + S*pi(-1) + eM;                //4.25
g0(13,rn) =1;
g0(13,mp_shock)=-1;
g0(13,pi)=-S;


% g = rhoG*g(-1) + eG;                            //4.26
g0(14,g)=1;
g1(14,g)=rhoG;
PSI(14,2)=sigma_G;
 
 %a = rhoA*a(-1) + eA;                            //4.27
g0(15,a)=1;
g1(15,a)=rhoA;
PSI(15,3)=sigma_A;

%Epi

g0(16,pi)=1;
g1(16,Epi)=1;
pi_matrix(16,1)=1;

%Ec

g0(17,c)=1;
g1(17,Ec)=1;
pi_matrix(17,2)=1;

%Erk

g0(18,rk)=1;
g1(18,Erk)=1;
pi_matrix(18,3)=1;


%mp shock
g0(19,mp_shock)=1;
g1(19,mp_shock)=rho;
PSI(19,1)=sigma_MP;

%call gensys
 
 [C,Constant,D_sqr,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,C,PSI,pi_matrix);
 D=D_sqr*D_sqr';

  %now add intercept term
 C_final=zeros(20,20);
 C_final(1:19,1:19)=C;
 C_final(20,20)=1;
 D_final=zeros(size(D,1)+1,size(D,2)+1);
 D_final(1:size(D,1),1:size(D,2))=D;
 
 
 
 
 A=zeros(3,20); %data must be in NOT annualized terms
 A(1,pi)=100; %inflation
 
 A(2,rn)=100;
 A(3,y)=100;% output gap

 B=zeros(3,3);



add_matrices=[];