function [f,c0]=Cem(z,m,q,ntrms)
% [f,c0]=Cem(z,m,q,ntrms)

% Modified Mathieu function Ce(z,m,q)
% z     - vector of argument values
% m     - integer function order
% q     - scalar parameter value
% ntrms - number of series terms used
% f     - vector of function values
% See 20.6.3 and 20.6.4 of 
% Abramowitz&Stegun.
if q==0
  f=cosh(m*z(:)); c0=[]; 
  if m==0, f=f/sqrt(2); end
  return;
end
if nargin<4, ntrms=100; end
z=z(:); nz=length(z); q2=sqrt(q);
z1=q2*exp(-z); z2=q2*exp(z);
if mod(m,2)==0 % even index   
  k=0:ntrms-1; sgn=cos(pi*k);    
  [a,c]=matue(q,1,1,ntrms);
  u=c(:,1+m/2); 
  c0=sum(u)*(sgn*u)/u(1)^2;
  f=(bes(k,z1).*bes(k,z2))*(c0*sgn'.*u);
else           % odd index
  k=0:ntrms-1; sgn=cos(pi*k);    
  [a,c]=matue(q,2,1,ntrms);
  u=c(:,1+fix(m/2)); 
  c0=sum(u)*(sgn.*(2*k+1))*u;
  c0=c0/(q2*u(1)^2);
  f=(bes(k,z1).*bes(k+1,z2)+...
     bes(k+1,z1).*bes(k,z2))*...
     (c0*sgn'.*u);
end

%=======================================

function f=bes(n,z)
% f=bes(n,z) calls function besselj
% for vector values of n and z.
% f(i,j)= besselj(n(j),z(i))
f= besselj(ones(size(z(:)))*n(:)',z(:)*ones(size(n(:)')));