function [A B C_final D_final add_matrices] = CK3( params,setup,data )
%code solves the model  the JME (2008) paper
%by Christoffel and Kuester. All data are detrended before estimation, hence there are no constants floating around.
% Most of this is taken from the dynare file in Wieland's macro model database. The MP rule was changed to depend only on
%one lag of inflation.


%unwrap parameters
kappa_p=params(1);
gamma_b=params(2);
epsilon = params(3);                % own price elasticity of demand for differentiated good
habit   = params(4);               % habit persistence
sig     = params(5);               % CRRA
vphi    = params(6);                 % inverse of Frisch elasticity


gamma       = params(7);  % Calvo wage stickiness  
wage_index  = params(8);    % wage indexation
vtheta      = params(9); % rate of separation


gamma_R     = params(10);	 % interest rate smoothing
gamma_Pi    = params(11);   % weight on inflation 
gamma_y     = params(12);   % weight on output

rho_eb      = params(13);   % persistence of preference shock
rho_emoney  = params(14);     % persistence of mon pol shock
rho_g       = params(15); % persistence of gov spending shock CP
rho_z       = params(16); % persistence of productivity shock CP

sig_innoeb =   params(17);
sig_monpol =   params(18); 
sig_innog  =   params(19);
sig_innoz  =   params(20);
%sig_u      =   params(21); %measurement error in inflation







%fixing other coefficients 
bet     = 0.9939; %time-discount factor
gamma_f=bet*(1-bet*gamma_b);
xi          = .5;  % weight on unemployment in matching function 
eta         = .5;  % bargaining power of workers 
alp         = .99; % labor elasticity of production 
%steady states

ybar    = 1;                    % steady state output
hbar    = 1/3;                  % steady state hours worked
Phi_y   = 0.00863;              % fraction of output lost due to fixed costs.
ubar    = 0.1;                  % steady state unemployment rate  
qbar    = 0.7;                  % probability of finding a worker in a given month
gbar    = 0.347026648444032;    % share of government spending in GDP (y=c+g)
b_wh    = 0.4;                  % replacement rate
Pibar   = 1;                    % zero inflation steady state


% ss for the given targets
nbar        = 1-ubar;   
mbar        = vtheta*nbar;
vbar        = mbar/qbar;
sbar        = mbar/ubar;
sigmam      = mbar*(ubar^xi*vbar^(1-xi))^(-1);
thetabar    = vbar/ubar;

mcbar       = (epsilon-1)/epsilon;
xLbar       = mcbar;
Rbar        = 1/bet;

zbar        = ybar/(nbar*hbar^(alp));

wbar        = xLbar*zbar*alp*hbar^(alp-1);
Phi         = Phi_y*ybar/nbar;

PsiLbar     = xLbar*zbar*hbar^(alp) - wbar*hbar - Phi;
PsiCbar     = (1-mcbar)*ybar;

Jbar		= 1/(1-bet*(1-vtheta))*PsiLbar;        
deltaFbar   = 1/(1-bet*(1-vtheta)*gamma)*wbar*hbar;

b           = b_wh*wbar*hbar;

mrsbar      = (Jbar*eta/(1-bet*(1-vtheta)*gamma)*alp/(alp-1)*hbar*wbar - (1-eta)*deltaFbar/(1-bet*(1-vtheta-sbar))*(wbar*hbar-b))/ (Jbar*hbar*eta/(1-bet*(1-vtheta)*gamma)*(-1)/(1-alp) - (1-eta)*deltaFbar*hbar/(1-bet*(1-vtheta-sbar))*1/(1+vphi)  );

deltaWbar   = 1/(1-bet*(1-vtheta)*gamma)*hbar*1/(1-alp)*(-alp*wbar - (-1)*mrsbar);
Deltabar    = 1/(1-bet*(1-vtheta-sbar))*(wbar*hbar - mrsbar*hbar/(1+vphi) - b); 

kappa       = qbar*bet*Jbar;
cbar        = ybar - kappa*vbar - Phi*nbar-gbar;

lambdabar   = (cbar*(1-habit))^(-sig);
kappaL      = mrsbar*lambdabar/hbar^vphi;



%indices for positions of variables

ct=1;
deltaFt=2;
deltaWt=3;
Deltastart=4;
ht=5;
Jstart=6;
lambdat=7;
mct=8;
mt=9;
nt=10;
Pit=11;
qt=12;
Rt=13;
st=14;
ut=15;
vt=16;
wstart=17;
wt=18;
xLt=19;
yt=20;
ebt=21;
emoneyt=22; 
gt=23;
zt=24; 
Epi=25;
Ewstar=26;
EdeltaF=27;
Elambda=28;
EdeltaW=29;
EJstar=30;
EDeltastar=31;
   
%matrices for gensys

g0=zeros(24+7,24+7);
g1=zeros(24+7,24+7);
c=zeros(24+7,1);
pi_matrix=zeros(24+7,7);
psi=zeros(24+7,4);


%model

%Monetary policy rule in the paper


g0(1,Rt)=1;
g0(1,yt)=-(1-gamma_R)*gamma_y/4;
g0(1,emoneyt)=-1;

g1(1,Pit)=(1-gamma_R)*gamma_Pi;
g1(1,Rt)=gamma_R;


% marginal utility of consumption
   

g0(2,lambdat)=1;
g0(2,ct)=sig/(1-habit);
g1(2,ct)=sig/(1-habit)*habit;

% New Keynesian Phillips curve (with inflation indexation)

g0(3,Pit)=1;
g0(3,Epi)=-gamma_f;
g0(3,mct)=-kappa_p;
g1(3,Pit)=gamma_b;

%marginal cost


g0(4,mct)=1;
g0(4,xLt)=-1;
     
	 
	 
% new matches

g0(5,mt)=1;
g0(5,ut)=-xi;
g0(5,vt)=-(1-xi);

%employment

 g0(6,nt)=1;
g1(6,nt)= (1-vtheta);
g1(6,mt)= mbar/nbar;


% unemployment 
                        
g0(7,nt)=1;
g0(7,ut)=ubar/(1-ubar);


% job-filling rate


g0(8,qt)=1;
g0(8,mt)=-1;
g0(8,vt)=1;

                            
% job-finding rate

g0(9,st)=1;
g0(9,mt)=-1;
g0(9,ut)=1;

%newly optimized wage (wage setting FOC)

       
g0(10,Jstart)=1;
g0(10,deltaWt)=1;
g0(10,Deltastart)=-1;
g0(10,deltaFt)=-1;


% hours FOC


g0(11,wt)=1;
g0(11,xLt)=-1;
g0(11,zt)=-1;
g0(11,ht)=-(alp-1);

% evolution of aggregate wage
 

g0(12,wt)=1;
g0(12,Pit)=gamma; 
g0(12,wstart)=-(1-gamma);
g1(12,wt)=gamma;
g1(12,Pit)=gamma*wage_index;


% deltaFt (-\partial surplus of firm/\partial wage)
	
g0(13,deltaFt)=1;	
g0(13,wstart)=(1-bet*(1-vtheta)*gamma)* (alp/(1-alp))+bet*(1-vtheta)*gamma*alp/(1-alp);
g0(13,xLt)=-(1-bet*(1-vtheta)*gamma)*1/(1-alp);
g0(13,zt)=-(1-bet*(1-vtheta)*gamma)*1/(1-alp);
g0(13,Pit)= bet*(1-vtheta)*gamma*alp/(1-alp)*wage_index;
g0(13,Ewstar)=-bet*(1-vtheta)*gamma*alp/(1-alp);
g0(13,Epi)=-bet*(1-vtheta)*gamma*alp/(1-alp);
g0(13,EdeltaF)=-bet*(1-vtheta)*gamma;
g0(13,Elambda)=-bet*(1-vtheta)*gamma;
g0(13,lambdat)=bet*(1-vtheta)*gamma;



% deltaWt ( \partial surplus of worker/\partial wage) 


g0(14,deltaWt)=deltaWbar;
g0(14,wstart)=-alp/(1-alp)*wbar*hbar*alp/(1-alp)+1/(1-alp)*mrsbar*hbar*(1+vphi)/(1-alp)-bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*( (alp/(1-alp))^2*wbar*hbar - (1+vphi)/(1-alp)^2*mrsbar*hbar );
g0(14,xLt)=alp/(1-alp)*wbar*hbar*1/(1-alp)-1/(1-alp)*mrsbar*hbar*(1+vphi)/(1-alp);
g0(14,zt)=alp/(1-alp)*wbar*hbar*1/(1-alp)-1/(1-alp)*mrsbar*hbar*(1+vphi)/(1-alp);
g0(14,lambdat)=1/(1-alp)*mrsbar*hbar+bet*(1-vtheta)*gamma*deltaWbar;
g0(14,Pit)=-bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*( (alp/(1-alp))^2*wbar*hbar - (1+vphi)/(1-alp)^2*mrsbar*hbar )*wage_index;
g0(14,Ewstar)=bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*( (alp/(1-alp))^2*wbar*hbar - (1+vphi)/(1-alp)^2*mrsbar*hbar );
g0(14,Epi)=bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*( (alp/(1-alp))^2*wbar*hbar - (1+vphi)/(1-alp)^2*mrsbar*hbar );
g0(14,Elambda)=-bet*(1-vtheta)*gamma*deltaWbar;
g0(14,EdeltaW)=-bet*(1-vtheta)*gamma*deltaWbar;


% value of firm that resets wage
             
					
g0(15,Jstart)=Jbar;
g0(15,wstart)=wbar*hbar+bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar;
g0(15,xLt)=	-wbar*hbar/alp;
g0(15,zt)=	-wbar*hbar/alp;	
g0(15,Ewstar)=-bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar;
g0(15,Epi)=-bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar;
g0(15,Pit)=bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar*wage_index;	
g0(15,Elambda)=	-bet*(1-vtheta)*Jbar;
g0(15,EJstar)=	-bet*(1-vtheta)*Jbar;
g0(15,lambdat)=	bet*(1-vtheta)*Jbar;



% surplus of worker that resets wage


g0(16,Deltastart)=Deltabar;
g0(16,wstart)= wbar*hbar*1/(1-alp)*alp-mrsbar*hbar/(1+vphi)*(1+vphi)/(1-alp)+bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar);
g0(16,xLt)= wbar*hbar*1/(1-alp)+mrsbar*hbar/(1+vphi)*(1+vphi)/(1-alp);
g0(16,zt)= wbar*hbar*1/(1-alp)+mrsbar*hbar/(1+vphi)*(1+vphi)/(1-alp); 
g0(16,lambdat)=-mrsbar*hbar/(1+vphi)+bet*(1-vtheta-sbar)*Deltabar;
g0(16,Pit)=bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar)*wage_index-bet*gamma*sbar/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar)*wage_index;
g0(16,Epi)=-bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar)+bet*gamma*sbar/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar);
g0(16,Ewstar)=-bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar)+bet*gamma*sbar/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar);
g0(16,wt)=-bet*gamma*sbar/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar);
g0(16,Elambda)=-bet*(1-vtheta-sbar)*Deltabar;
g0(16,EDeltastar)=-bet*(1-vtheta-sbar)*Deltabar;
g0(16,st)=bet*Deltabar*sbar;


% vacancy posting condition

   g0(17,qt)=-kappa/qbar;
   g0(17,wt)=bet*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar;
   g0(17,Pit)=bet*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar*wage_index;
   g0(17,Ewstar)=-bet*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar;
    g0(17,Epi)=-bet*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar;
	g0(17,Elambda)=-bet*Jbar;
	g0(17,EJstar)=-bet*Jbar;
	g0(17,lambdat)=bet*Jbar;
	
	
	% resource constraint
 
 
 g0(18,yt)=ybar;
 g0(18,ct)=-cbar;
 g0(18,gt)=-gbar;
 g0(18,vt)=-vbar*kappa;
 g0(18,nt)=-nbar*Phi;
 
 
% production
   
 g0(19,yt)=1;
g0(19,nt)=-1; 
g0(19,zt)=-1;
g0(19,ht)=-alp;


% ****************** shocks ************************
% shock to discount factor

 
g0(20,ebt)=1;
g1(20,ebt)=rho_eb;
psi(20,1)=sig_innoeb; 


% government spending

g0(21,gt)=1;
g1(21,gt)=rho_g;
psi(21,2)=sig_innog;    
              
% monetary policy shock
g0(22,emoneyt)=1;
g1(22,emoneyt)=rho_emoney;
psi(22,3)=sig_monpol;   

 
%productivity shock

g0(23,zt)=1;
g1(23,zt)=rho_z;
psi(23,4)=sig_innoz;


%consumption Euler equation
    

g0(24,lambdat)=1;
g0(24,Elambda)=-1;
g0(24,Rt)=-1;
g0(24,ebt)=-1;
g0(24,Epi)=1;


%definition of expectations

%pi

g0(25,Pit)=1;
g1(25,Epi)=1;
pi_matrix(25,1)=1;

%wstar

g0(26,wstart)=1;
g1(26,Ewstar)=1;
pi_matrix(26,2)=1;

%deltaF

g0(27,deltaFt)=1;
g1(27,EdeltaF)=1;
pi_matrix(27,3)=1;

%lambda

g0(28,lambdat)=1;
g1(28,Elambda)=1;
pi_matrix(28,4)=1;

%deltaW

g0(29,deltaWt)=1;
g1(29,EdeltaW)=1;
pi_matrix(29,5)=1;


%Jstar

g0(30,Jstart)=1;
g1(30,EJstar)=1;
pi_matrix(30,6)=1;

%Deltastar

g0(31,Deltastart)=1;
g1(31,EDeltastar)=1;
pi_matrix(31,7)=1;


%call gensys
 
 [C,Constant,D_sqr,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi_matrix);
 D=D_sqr*D_sqr';

 %now add intercept term
 C_final=zeros(32,32);
 C_final(1:31,1:31)=C;
 C_final(32,32)=1;
 D_final=zeros(size(D,1)+1,size(D,2)+1);
 D_final(1:size(D,1),1:size(D,2))=D;
 
 
 
 
 A=zeros(3,32);
 A(1,Pit)=100; %inflation
 
 A(2,Rt)=100;
 A(3,yt)=100;% output gap
 A(4,wt)=100;%wage
 %A(5,ut)=100;%unemployment
 B=zeros(4,4);
 %B(5,5)=sig_u^2;



add_matrices=[];