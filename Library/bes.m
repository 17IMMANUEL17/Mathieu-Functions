function f=bes(n,z)
% f=bes(n,z) calls function besselj
% for vector values of n and z.
% f(i,j)= besselj(n(j),z(i))
f= besselj(ones(size(z(:)))*n(:)',z(:)*ones(size(n(:)')));