function f=Ms1(z,m,q,ntrms)
% f=Ms1(z,m,q,ntrms)

% Modified Mathieu function Ms1(z,m,q)
% defined by 20.6.9 & 20.6.10 of
% 'Handbook of Mathematical Functions'
% by Milton Abramowitz and Irene Stegun
% z     - vector of argument values
% m     - integer function order
% q     - scalar parameter value
% ntrms - number of series terms used
% f     - vector of function values

%       by H.W. 12/01/04

if nargin<4, ntrms=50; end
z=z(:); z1=sqrt(q)*exp(-z); 
z2=sqrt(q)*exp(z); 
if mod(m,2)==0 % even order
  k=1:ntrms; [a,c]=matue(q,1,2,ntrms);
  u=c(:,m/2); r=m/2;  
  [ubig,s]=max(abs(u)); p0=(-1)^r/u(s);
  f=(bes(k-s,z1).*bes(k+s,z2)-...
     bes(k+s,z1).*bes(k-s,z2))*...
     (p0*cos(k*pi)'.*u);
else           % odd order
  k=0:ntrms-1; [a,c,vsave]=matue(q,2,2,ntrms);
  u=c(:,1+fix(m/2)); [ubig,s]=max(abs(u));
  r=(m-1)/2; p0=(-1)^r/u(s); s=s-1;
  f=(bes(k-s,z1).*bes(k+s+1,z2)-...
     bes(k+s+1,z1).*bes(k-s,z2))*...
     (p0*cos(k*pi)'.*u);
end
end

%========================================
% 
% function f=bes(n,z)
% % f=bes(n,z) calls function besselj
% % for vector values of n and z.
% % f(i,j)= besselj(n(j),z(i))
% f= besselj(ones(size(z(:)))*n(:)',z(:)*ones(size(n(:)')));