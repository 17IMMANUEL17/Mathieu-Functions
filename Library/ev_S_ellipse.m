%% Eigenvalues of S on the ellipse
% See SPECTRAL DECOMPOSITIONS OF BOUNDARY INTEGRAL OPERATORS
clear all
close all
clc
nMax = 110; 
a = 1; % ellipse param 1
mu = 1.2; % ellipse param 2
a_x = a * cosh(mu);
b_x = a * sinh(mu);
k = 2;
qP = 1/4 * (k * a)^2; % Parameter ellipse
load('sv_ordered_mu1-2_k4_a1_401.mat')
% Spectral decompositions and nonnormality of boundary integraloperators
% in acoustic scattering - extended version,T. BETCKE†, J. PHILLIPS‡, 
% E. A. SPENCE§: equation before 4.3 
Mc3 = @(n,mu,q) Mc1(mu,n,q,nMax) + 1i * Mc2(mu,n,q,nMax);
Ms3 = @(n,mu,q) Ms1(mu,n,q,nMax) + 1i * Ms2(mu,n,q,nMax);
% Eigenvalue decomposition S: equations after 4.3
ldaC = @(n,mu,q) ldaC_S(n,mu,q,nMax);
ldaS = @(n,mu,q) ldaS_S(n,mu,q,nMax);
%as= @(n,mi,q) n*exp(n*mi)*Mc2(mi,n,q)/factorial(n)*2^(-n+1)*(q^(n/2));
a=1;
k=4;
q=1/4*(k*a)^2;

% for n=1:200
% b(n)=ldaC_S(n-1,1.1,1,nMax);
% s(n)=ldaS_S(n,1.1,1,nMax);
% end
% b(201)=ldaC_S(200,1.1,1,nMax);
% s=fliplr(s);
% se=cat(2,s,b);
% se=abs(se);
% figure
%  plot((1:1:401),se,'xb');
%  hold on
%  c=singularValuesS/1.5067;
%  plot((1:1:401),c,'or')
%  legend('analytic eigenvalues','numerical eigenvalues')
%  figure
%  err_rel= abs(se'-c)./abs(c);
%  plot((1:401),err_rel)
 

 for n=1:200
b(n)=ldaC_S(n-1,mu,q,nMax);
s(n)=ldaS_S(n,mu,q,nMax);
end
b(201)=ldaC_S(200,mu,q,nMax);
s=fliplr(s);
se=cat(2,s,b);
se=abs(se);
 plot((1:1:401),se,'b');
 hold on
 c=sv_ordered/radius(1,1.2)%1.6635;
 plot((1:1:401),c','r')
 legend('analytic eigenvalues','numerical eigenvalues')
 figure
 err_rel= abs(se'-c)./abs(c);
 semilogy((1:401),err_rel)
 legend('err rel')
%    

% %% Get Analytic solution (not working)
% for nn = 1 : 20
%     ldaCVec(nn) = ldaC(nn,1,q) ;
%     ldaSVec(nn) = ldaS(nn,1,q) ;
% end
% 
% theta = linspace(-pi,pi,100);
% Jac = sqrt((a_x * sin(theta)).^2 + (b_x * cos(theta)).^2);
% 
% phi_i = 3*pi/2;
% j = zeros(1,length(theta));
% m = 0;
% nEv = 2;
% for n = -nEv:nEv
%     m = m+1;
%     Beta_C_n = pw_TM_EFIE_ellipse_coeff_C(phi_i, n,mu,q,nMax);
%     Beta_S_n = pw_TM_EFIE_ellipse_coeff_S(phi_i, n,mu,q,nMax);
%     if n ~= 0        
%         alpha_C_n= ldaC(abs(n),mu,q);
%         alpha_S_n= ldaS(abs(n),mu,q);
%         j=j+ ce(theta,abs(n),q).'*Beta_C_n/alpha_C_n./Jac;      
%         j=j+ se(theta,abs(n),q).'*Beta_S_n/alpha_S_n./Jac;   
%     else
%         alpha_C_n= ldaC(n,mu,q);       
%         j=j+ ce(theta,n,q).' *Beta_C_n/alpha_C_n./Jac;        
%     end                       
% end