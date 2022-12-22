%% Eigenvalues of K on the ellipse
% See SPECTRAL DECOMPOSITIONS OF BOUNDARY INTEGRAL OPERATORS
nMax = 50; % maximum order
a = 1; % ellipse param 1
mu = 3; % ellipse param 2
a_x = a * cosh(mu);
b_x = a * sinh(mu);
k = 2;
qP = 1/4 * (k * a)^2; % Parameter ellipse

%  equation before 4.3 
Mc3 = @(n,mu,q) Mc1(mu,n,q,nMax) + 1i * Mc2(mu,n,q,nMax);
Ms3 = @(n,mu,q) Ms1(mu,n,q,nMax) + 1i * Ms2(mu,n,q,nMax);
% Eigenvalue decomposition K
ldaC = @(n,mu,q) 1i*pi/2 * Mc1p(mu,n,q,nMax) * Mc3(n,mu,q);
ldaS = @(n,mu,q) 1i*pi/2 * Ms1p(mu,n,q,nMax) * Ms3(n,mu,q);
