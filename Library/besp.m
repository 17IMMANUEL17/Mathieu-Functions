function fp=besp(n,z)
% fp=besp(n,z) computes the derivative of 
% besselj for vector values of n and z
n=n(:)'; z=z(:);
fp=(bes(n-1,z)-bes(n+1,z))/2;