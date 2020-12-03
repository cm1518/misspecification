function [A B C_final D_final add_matrices] = HS(para,setup,data)
% Solve the Small-Scale DSGE model in the textbook 
% INPUT
% para:  structural parameters
%
% OUTPUT
% para:  structural parameters
%
% DATE: 2/13/2016
% Written by Minsu Chang


%=========================================================================
%                     Paramaters
%=========================================================================


kappa = para(1);
rho_z = para(2);
tau = para(3);
psi1 = para(4);
psi2 = para(5);
rA = para(6);
piA = para(7);
gammaQ = para(8);
rho_R = para(9);
rho_g = para(10);

sigma_R = para(11);
sigma_g = para(12);
sigma_z = para(13);

bet = 1/(1+rA/400);


%=========================================================================
%                        SOLVE DSGE MODEL
%=========================================================================


retcode = 0;

valid   = 1;

% /* Variable indices */

%=========================================================================
%                       DEFINE OBJECTS
%=========================================================================

% Equation indices

eq_1    = 1;  %** (2.1) on \hat{y}(t) **/
eq_2    = 2;  %** (2.1) on \hat{pi}(t) **/
eq_3    = 3;  %** (2.1) on \hat{R}(t) **/
eq_4    = 4;  %** \hat{y}(t-1) **/
eq_5   = 5;  %** \hat{g} process **/
eq_6    = 6;  %** \hat{z} process **/
eq_7 = 7;  %** \hat{y} expectational error **/
eq_8   = 8;  %** \hat{pi} expectational error **/

% Variable indices 

y_t   = 1;
pi_t   = 2;
R_t   = 3;
y1_t  = 4;
g_t   = 5;
z_t   = 6;
Ey_t1   = 7;
Epi_t1  = 8;

%Expectation error indices (eta) 

ey_sh  = 1;
epi_sh  = 2;
       
%Shock indices (eps)

z_sh = 1;
g_sh = 2;
R_sh = 3;

%SUMMARY

neq  = 8;
neta = 2;
neps = 3;

% /** initialize matrices **/

GAM0 = zeros(neq,neq);
GAM1 = zeros(neq,neq);
   C = zeros(neq,1);        
 PSI = zeros(neq,neps);
 PPI = zeros(neq,neta);



%=========================================================================
%                EQUILIBRIUM CONDITIONS: CANONICAL SYSTEM
%=========================================================================

%=========================================================================
%         1. 
%=========================================================================

GAM0(eq_1,y_t)   =  1;
GAM0(eq_1,R_t)   =  1/tau;
GAM0(eq_1,g_t)  = -(1-rho_g);
GAM0(eq_1,z_t)  = -rho_z/tau;
GAM0(eq_1, Ey_t1) = -1;
GAM0(eq_1, Epi_t1) = -1/tau;

%=========================================================================
%         2. 
%=========================================================================

GAM0(eq_2,y_t)   =  -kappa;
GAM0(eq_2,pi_t)   = 1;
GAM0(eq_2,g_t)   =  kappa;
GAM0(eq_2, Epi_t1) = -bet;

%=========================================================================
%         3. 
%=========================================================================

GAM0(eq_3,y_t)   = -(1-rho_R)*psi2;
GAM0(eq_3,pi_t)  = -(1-rho_R)*psi1;
GAM0(eq_3,R_t)  = 1;
GAM0(eq_3,g_t) = (1-rho_R)*psi2;
GAM1(eq_3,R_t) = rho_R;
PSI(eq_3,R_sh) = sigma_R;
  
%=========================================================================
%         4. 
%=========================================================================

GAM0(eq_4,y1_t)   = 1;
GAM1(eq_4,y_t) = 1;

%=========================================================================
%         5. 
%=========================================================================

GAM0(eq_5,g_t)   = 1;
GAM1(eq_5, g_t) = rho_g;
PSI(eq_5, g_sh) = sigma_g;

%=========================================================================
%         6. 
%=========================================================================

GAM0(eq_6,z_t)   = 1;
GAM1(eq_6, z_t) = rho_z;
PSI(eq_6, z_sh) = sigma_z;

%=========================================================================
%         7. 
%=========================================================================

GAM0(eq_7, y_t) = 1;
GAM1(eq_7, Ey_t1) = 1;
PPI(eq_7, ey_sh) = 1;

%=========================================================================
%         8. 
%=========================================================================

GAM0(eq_8, pi_t) = 1;
GAM1(eq_8, Epi_t1) = 1;
PPI(eq_8, epi_sh) = 1;

      

%=========================================================================
%           QZ(generalized Schur) decomposition by GENSYS
%=========================================================================

[C,Constant,D_sqr,fmat,fwt,ywt,gev,eu,loose] = gensys(GAM0,GAM1,C,PSI,PPI,1+1E-8);
 D=D_sqr*D_sqr';

 %now add intercept term
 C_final=zeros(size(C,1)+1,size(C,2)+1);
 C_final(1:size(C,1),1:size(C,2))=C;
 C_final(end,end)=1;
C_final(1:end-1,end)=Constant;


 D_final=zeros(size(D,1)+1,size(D,2)+1);
 D_final(1:size(D,1),1:size(D,2))=D;
 

 
 
 A=zeros(3,9);
 A(1,2)=1; %inflation
 A(1,9)=piA/4;
 A(2,3)=1;
A(2,9)=piA/4+rA/4+gammaQ; %nominal interest rate
 A(3,1)=1;% GDP growth
 A(3,4)=-1;% GDP growth
 A(3,6)=1;% GDP growth
 A(3,9)=gammaQ;
 B=zeros(3,3);



add_matrices=[];



end
