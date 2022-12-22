%% Eigenvalues of T on the ellipse
% See SPECTRAL DECOMPOSITIONS OF BOUNDARY INTEGRAL OPERATORS
nMax = 50; % maximum order
a = 1; % ellipse param 1
mu = 3; % ellipse param 2
a_x = a * cosh(mu);
b_x = a * sinh(mu);
k = 2;
qP = 1/4 * (k * a)^2; % Parameter ellipse

% Eigenvalue decomposition T
ldaC = @(n,mu,q) 1i*pi/2 * Mc1(mu,n,q,nMax) * (Mc1p(mu,n,q,nMax) + 1i * Mc2p(mu,n,q,nMax));
ldaS = @(n,mu,q) 1i*pi/2 * Ms1(mu,n,q,nMax) * (Ms1p(mu,n,q,nMax) + 1i * Ms2p(mu,n,q,nMax));
ldaS(1,1,10)